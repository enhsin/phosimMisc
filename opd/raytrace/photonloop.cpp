///
/// @package phosim
/// @file photonloop.cpp
/// @brief main photon loop
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "event.h"
#include "grating.h"
#include "counter.h"

int Image::photonLoop () {

    Vector angle, largeAngle, position, positionDiffraction;
    Vector positionPrevious, anglePrevious, normal;
    double transmission, reflection, distance;

    long long ray;
    long long photonSource;
    long sourceCounter;

    long long detRay;
    long long backSourceOver;
    double sourceSaturation;
    double backSigma = 0.0;
    long preGhost = 0;
    long sourceOver;
    long newSurf = 0;
    long oldSurf;
    long waveIndex;
    long waveSurfaceIndex;
    long straylightcurrent;
    long surfaceLimit;

    int miss;
    int initFlag;
    char tempstring[4096];

    double cx, cy, cz, r0;
    int writeEvent=eventfile;
    long targetRay = 80*OPD_SCREEN_SIZE + 15;

    EventFile* pEventLogging;
    Clog counterLog;
    Tlog throughputLog;

    Grating* pGrating;

    //    OPTIMIZATION INITIALIZATION
    photon.maxcounter = 0;
    surfaceLimit = natmospherefile*2 + nsurf*2 + 2;
    dynamicTransmission = (double*)malloc(surfaceLimit*901*sizeof(double));
    for (int i = 0; i < 901; i++) {
        for (int j = 0; j < surfaceLimit; j++) {
            dynamicTransmission[i*surfaceLimit + j] = -1.0;
        }
    }
    if (checkpointcount != 0) readCheckpoint(checkpointcount);

    //    LOGGING INITIALIZATION
    if (eventfile) {
        pEventLogging = new EventFile((int)(nphot*100), outputdir, eventFitsFileName);
    } else {
        pEventLogging = NULL;
    }
    pGrating = new Grating();
    if (throughputfile) initThroughput(&throughputLog, nsurf);
    counterInit(&counterLog);

    //    SOURCE TYPE LOOP
    for (int sourceType = 0; sourceType < 243; sourceType++) {
        sourceCounter = 0;

        if (static_cast<int>(round(sourceType*checkpointtotal/242)) == checkpointcount) {

            //    SOURCE LOOP
            for (long source = 0; source < nsource; source++) {

                if ((sourceType >= 0 && sourceType <= 49 && sources.type[source] == 0 && (source%50) == sourceType) ||
                    (sourceType >= 50 && sourceType <= 99 && sources.type[source] == 1 && (source%50) == (sourceType - 50)) ||
                    (sourceType >= 100 && sourceType <= 149 && sources.type[source] == 2 && (source%50) == (sourceType - 100)) ||
                    (sourceType >= 150 && sourceType <= 199 && sources.type[source] == 3 && (source%50) == (sourceType - 150)) ||
                    (sourceType == 200 && sources.type[source] == 4 && sources.mag[source] > 40.5) ||
                    (sourceType >= 201 && sourceType <= 241 && sources.type[source] == 4 && sources.mag[source] > (240.5 - sourceType)
                     && sources.mag[source] <= (241.5 - sourceType))
                    || (sourceType == 242 && sources.type[source] == 4 && sources.mag[source] <= -0.5)) {

                    //

                    //    SETUP PHOTONS
                    photonSource = static_cast<long long>(nphot*sources.norm[source]/totalnorm);
                    if (nsource == 1) photonSource = nphot;
                    photonSource = poisson(photonSource);
                    if (telconfig != 0 && sources.type[source] != 0) photonSource = 0;
                    if (opdfile) photonSource = OPD_SCREEN_SIZE*OPD_SCREEN_SIZE;

                    ray = 0;
                    sourceOver = 1;
                    photon.sourceOver_m = 1;
                    sourceSaturation = 1.0;
                    photon.sourceSaturationRadius = 0.0;
                    obstruction.pupil = 0;
                    photon.prtime = -1.0;
                    if (sources.type[source] < 4) {
                        double rate = (static_cast<double>(counterLog.accepted) + 10)/(static_cast<double>(counterLog.totalPhoton) + 1000);
                        if (rate < 0.01) rate = 0.01;
                        if (backGamma < 1.0) backSourceOver = (long long)(backAlpha*sqrt(rate*photonSource));
                        else backSourceOver = (long long)(backAlpha/backGamma*sqrt(rate*photonSource));
                        if (backAlpha <= 0.0) backSourceOver = static_cast<long long>(backGamma);
                        if (backSourceOver < 1) backSourceOver = 1;
                        if (backSourceOver > photonSource) backSourceOver = photonSource;
                    } else {
                        backSourceOver = 1;
                    }
                    if (sources.type[source] >= 4) {
                        backSigma = 0.0;
                    } else {
                        backSigma = backBeta*backRadius/3600.0*platescale/1000.0;
                    }
                    if (sources.mag[source] < straylightcut && straylight == 1) {
                        straylightcurrent = 1;
                    } else {
                        straylightcurrent = 0;
                    }

                    //   MAIN PHOTON LOOP
                photonloop:
                    while (ray < photonSource) {
                        if (writeEvent && (ray==0 || ray==targetRay)) eventfile=1; else eventfile=0;
                        sourceOver = 1;
                        photon.sourceOver_m = 1;

                        //   Get Wavelength and Time
                        miss = getWavelengthTime(&(photon.wavelength), &(photon.time), source);
                        photon.absoluteTime = photon.time - exptime/2.0 + timeoffset;
                        // Get Angle
                        getAngle(&angle, photon.time, source);
                        if (eventfile) {
                            pEventLogging->logPhoton(angle.x, angle.y, photon.wavelength, 0);
                            pEventLogging->logPhoton(photon.time, 0, 0, 1);
                        }

                        //   Saturation Optimization
                        if (saturation) {
                            if (sourceSaturation > 1.0 && photon.sourceSaturationRadius > 0.0 && sources.type[source] >= 4) {
                                sourceOver = round(sourceSaturation);
                                if (sourceOver > well_depth) sourceOver = well_depth;
                                if (sourceOver < 1) sourceOver = 1;
                                photon.sourceOver_m = round(((sourceSaturation - np)/(1 - np)));
                                if (photon.sourceOver_m < 1) {
                                    photon.sourceOver_m = 1;
                                    sourceOver = 1;
                                }
                            }
                        }
                        if (sources.type[source] < 4 && backGamma > 1.0) {
                            sourceOver = static_cast<long long>(backGamma);
                            photon.sourceOver_m = sourceOver;
                            if (backSourceOver*sourceOver > photonSource) {
                                sourceOver = static_cast<long long>(photonSource/backSourceOver);
                                photon.sourceOver_m = static_cast<long long>(photonSource/backSourceOver);
                            }
                            if (sourceOver < 1) {
                                sourceOver = 1;
                                photon.sourceOver_m = 1;
                            }
                        }


                        if (miss) {
                            countBad(&counterLog, sourceOver*backSourceOver, &ray);
                            goto photonloop;
                        }
                        waveIndex = static_cast<long>(photon.wavelength*1000 - 300);
                        waveSurfaceIndex = waveIndex*surfaceLimit;
                        initFlag = 0;
                        if (throughputfile) addThroughput(&throughputLog, -1, waveIndex, sourceOver*backSourceOver);

                        // Dynamic Transmission Optimization
                        long kmax;
                        if (sources.type[source] < 4 && backGamma > 1.0) {
                            kmax = 1;
                        } else {
                            kmax = sourceOver;
                        }
                        for (long k = 0; k < kmax; k++) {

                            if ((k == 0) || (k > 0 && straylightcurrent == 1)) {
                                long lastSurface = -1;
                                miss = dynamicTransmissionOptimization(k, &lastSurface, &preGhost, waveSurfaceIndex, straylightcurrent);
                                if (miss == 1) {
                                    if (throughputfile && lastSurface >= 0) {
                                        for (long j = 0; j <= lastSurface; j++) {
                                            addThroughput(&throughputLog, j, waveIndex, sourceOver*backSourceOver);
                                        }
                                    }
                                    countBad_dt(&counterLog, sourceOver*backSourceOver, &ray);
                                    goto photonloop;
                                }
                            }

                        redodiff:;
                            vectorInit(&largeAngle);

                            //  Sample Pupil
                            miss = samplePupil(&positionDiffraction, ray);
                            if (miss) {
                                countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                goto photonloop;
                            }
                            miss = samplePupil(&position, ray);
                            if (miss) {
                                countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                goto photonloop;
                            }
                            if (finiteDistance != 0.0) {
                                double finiteDistanceR = sqrt(pow(finiteDistance, 2) + pow(position.x, 2) +
                                                              pow(position.y, 2));
                                angle.x += position.x/finiteDistanceR;
                                angle.y += position.y/finiteDistanceR;
                                angle.z = smallAnglePupilNormalize(angle.x, angle.y);
                            }
                            if (diffraction_on == 5) vectorCopy(position, &positionDiffraction);
                            photon.xp = position.x;
                            photon.yp = position.y;
                            if (initFlag == 0) {
                                photon.shiftedAngle = spiderangle + photon.time*rotationrate*ARCSEC;
                                photon.wavelengthFactor = pow(photon.wavelength, -0.2)/screen.wavelengthfactor_nom;
                                vectorInit(&largeAngle);
                                initFlag = 1;
                            }

                            //  Diffraction
                            if (diffraction_on >= 1) {
                                miss = diffraction(&positionDiffraction, angle, &largeAngle);
                                if (miss) {
                                    if (k > 0) goto redodiff;
                                    countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                    goto photonloop;
                                }
                            }

                            //   Large Angle Scattering
                            largeAngleScattering(&largeAngle);

                            //   Second Kick
                            if (diffraction_on == 1 && sources.type[source] !=  0) secondKick(&largeAngle);

                            // Saturation Optimization Calculation
                            if (photon.sourceSaturationRadius > 0.0) {
                                if (modulus(&largeAngle)/DEGREE*platescale/pixsize
                                    > photon.sourceSaturationRadius || preGhost>= 2) {
                                    photon.sourceOver_m = 1;
                                    sourceSaturation -= ((1.0 - np)/np + 1.0)*0.001;
                                    if (sourceSaturation<1) sourceSaturation = 1;
                                    break;
                                }
                            }

                        }
                        if (photon.sourceSaturationRadius > 0.0) sourceSaturation += 0.001;
                        if (np <= 0.0) sourceSaturation = 1;
                        photon.counter = -1;
                        if (opdfile) photon.op = 0.0;

                        // Get Delta Angle
                        getDeltaAngle(&angle, &position, source);
                        photon.airRefraction = airIndexRefraction();
                        photon.ncurr = 1.0 + photon.airRefraction/1e6;

                        // ATMOSPHERE

                        // Atmospheric Dispersion
                        if (sources.type[source] !=  0) atmosphericDispersion(&angle);

                        // Loop over Atmosphere Layers
                        for (int layer = -1; layer < natmospherefile; layer++) {

                            // Atmosphere Propagate
                            atmospherePropagate(&position, angle, layer, diffraction_on);

                            if (sources.type[source] != 0) {
                                if (layer >= 0) {

                                    // Atmosphere Intercept
                                    atmosphereIntercept(&position, layer);

                                    // Atmosphere Refraction
                                    atmosphereRefraction(&angle, layer, diffraction_on);

                                    // Clouds
                                    transmission = cloudOpacity(layer);
                                    if (transmissionCheck(transmission, 1 + layer*2, waveSurfaceIndex)) {
                                        countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                        goto photonloop;
                                    }
                                }

                                // Atmosphere Opacity
                                transmission = atmosphereOpacity(angle, layer);
                                if (transmissionCheck(transmission, 2 + layer*2, waveSurfaceIndex)) {
                                    countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                    goto photonloop;
                                }

                                if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, layer + 100);

                            }
                        }

                        // Atmosphere Diffraction
                        if (diffraction_on == 2 && sources.type[source] != 0) atmosphereDiffraction(&angle);

                        // Dome Seeing
                        if (domeseeing > 0.0 || toypsf > 0.0) domeSeeing(&angle);

                        // Tracking
                        if (tracking_on) tracking(&angle, photon.time);

                        // Large Angle
                        angle.x += largeAngle.x;
                        angle.y += largeAngle.y;
                        angle.z = smallAnglePupilNormalize(angle.x, angle.y);

                        if (telescope_on == 0) {
                            newSurf = nsurf - 2;
                            if (nmirror % 2 == 0) {
                                propagate(&position, angle, ((surface.height[nsurf - 1] + platescale/DEGREE/1000) - position.z)/angle.z);
                                angle.x -= position.x/(platescale/DEGREE/1000);
                                angle.y -= position.y/(platescale/DEGREE/1000);
                                angle.z = -1.0;
                                normalize(&angle);
                            } else {
                                propagate(&position, angle, ((surface.height[nsurf - 1] - platescale/DEGREE/1000) - position.z)/angle.z);
                                angle.x -= position.x/(platescale/DEGREE/1000);
                                angle.y -= position.y/(platescale/DEGREE/1000);
                                angle.z = 1.0;
                                normalize(&angle);
                            }
                            if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, 200);
                        } else {
                            newSurf = -1;
                        }

                        // OPTICS AND DETECTOR
                        photon.direction = 1;
                        photon.ghostFlag = 0;
                    surfaceloop: while (1) {
                            oldSurf = newSurf;
                            if (photon.direction == 1) {
                                newSurf++;
                            } else {
                                newSurf--;
                            }

                        redostraylight:;

                            // Find intercept
                            if (newSurf>= 0 && newSurf<nsurf) {

                                transform(&position, &angle, newSurf);
                                if (surface.surfacetype[newSurf] == DETECTOR) {
                                    transform(&position, &angle, newSurf + 1);
                                }
                                miss = findSurface(angle, position, &distance, newSurf);

                            } else {
                                miss = 1;
                            }

                            //   Missed surface or wrong direction
                            //if ( eventfile ) printf("A: surf %d miss %d %f\n",newSurf,miss,distance);
                            //if ( eventfile ) printf("%g %g %g %g %g %g\n",position.x,position.y,position.z,angle.x,angle.y,angle.z);
                            if (miss == 1 || (distance < 0)) {
                                if (straylightcurrent == 0) {
                                    countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                    goto photonloop;
                                } else {
                                    if (chooseSurface(&newSurf, &oldSurf) == 1) {
                                        countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                        goto photonloop;
                                    } else {
                                        goto redostraylight;
                                    }
                                }
                            }

                            propagate(&position, angle, distance);


                            if (opdfile) {
                                if (newSurf == 0) {
                                    // assuming the entrance pupil at z=0. we starts calculating the optical path
                                    // from z=20m plane (just to be consistent with Bo's zemax calculation)
                                    photon.opdx = -position.z*angle.x/angle.z + position.x;
                                    photon.opdy = -position.z*angle.y/angle.z + position.y;
                                    distance = position.x*angle.x + position.y*angle.y + (position.z-20000)*angle.z;
                                }
                                photon.op -= distance*photon.ncurr;
                            }
                            if ( ray <= 1 ) printf("%d %g %g %g %g %g %g %g\n",newSurf,position.x,position.y,position.z,angle.x,angle.y,angle.z,photon.op);

                            if (throughputfile == 1 && photon.direction == 1) {
                                addThroughput(&throughputLog, newSurf, waveIndex, sourceOver*backSourceOver);
                            }

                            if (photon.direction == -1) photon.ghostFlag = 1;

                            //   DERIVATIVES
                            interceptDerivatives(&normal, position, newSurf);

                            //   CONTAMINATION
                            if (surface.surfacetype[newSurf] !=  DETECTOR && contaminationmode==1) {
                                miss = contaminationSurfaceCheck(position, &angle, newSurf);
                                if (miss) {
                                    countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                    goto photonloop;
                                }
                            }

                            //   SURFACE COATINGS
                            transmission = surfaceCoating(photon.wavelength, angle, normal, newSurf, &reflection);
                            //if ( ray == 0 ) printf("B: surf %d\n",newSurf);
                            //if ( ray == 0 ) printf("%g %g %g %g %g %g\n",position.x,position.y,position.z,angle.x,angle.y,angle.z);

                            if (transmissionCheck(transmission, natmospherefile*2 + 1 + newSurf*2, waveSurfaceIndex)) {
                                if (straylightcurrent == 1 && ghost[newSurf] == 0) {
                                    if (transmissionCheck(reflection + transmission, natmospherefile*2 + 1 + newSurf*2 + 1, waveSurfaceIndex)) {
                                        countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                        goto photonloop;
                                    } else {
                                        photon.direction = -photon.direction;
                                        reflect(&angle, normal);
                                        transformInverse(&position, &angle, newSurf);
                                        if (surface.surfacetype[newSurf] == DETECTOR) {
                                            transformInverse(&position, &angle, newSurf + 1);
                                        }
                                        goto surfaceloop;
                                    }
                                } else {
                                    countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                   goto photonloop;
                                }
                            }
                            //if ( ray == 0 ) printf("C: surf %d\n",newSurf);
                            //if ( ray == 0 ) printf("%g %g %g %g %g %g\n",position.x,position.y,position.z,angle.x,angle.y,angle.z);


                            //   INTERACTIONS
                            if (surface.surfacetype[newSurf] == MIRROR) {

                                //   MIRROR
                                reflect(&angle, normal);
                                //if ( ray == 0 ) printf("C1 %g %g %g %g %g %g\n",position.x,position.y,position.z,angle.x,angle.y,angle.z);
                                transformInverse(&position, &angle, newSurf);
                                if ( ray == targetRay && newSurf == 0) printf("C2 %.10f %.10f %.10f %g %g %g\n",position.x,position.y,position.z,angle.x,angle.y,angle.z);
                                if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);
                                if (eventfile && opdfile) pEventLogging->logPhoton(distance, photon.ncurr, photon.op, newSurf + 400);

                            } else if (surface.surfacetype[newSurf] == LENS || surface.surfacetype[newSurf] == FILTER) {

                                //   LENS/FILTER
                                newRefractionIndex(newSurf);
                                refract(&angle, normal, photon.nprev, photon.ncurr);
                                transformInverse(&position, &angle, newSurf);
                                if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);
                                if (eventfile && opdfile) pEventLogging->logPhoton(distance, photon.nprev, photon.op, newSurf + 400);

                            } else if (surface.surfacetype[newSurf] == GRATING) {

                                //   GRATING
                                double wavelengthNm = photon.wavelength*1000.0;
                                Vector angleOut;
                                pGrating->diffract(angle.x, angle.y, angle.z, normal.x, normal.y, normal.z,
                                                   angleOut.x, angleOut.y, angleOut.z, wavelengthNm);
                                vectorCopy(angleOut, &angle);
                                transformInverse(&position, &angle, newSurf);
                                if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);

                            } else if (surface.surfacetype[newSurf] == DETECTOR) {

                                if (eventfile || opdfile) {
                                    transformInverse(&position, &angle, newSurf + 1);
                                    transformInverse(&position, &angle, newSurf);
                                    if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);
                                    if (eventfile && opdfile) pEventLogging->logPhoton(distance, photon.ncurr, photon.op, newSurf + 400);
                                    if (opdfile) {
                                        //find exit pupil
                                       //if ( ray == 0 ) {
                                       //    double R0 = surface.outerRadius[newSurf + 2];
                                       //    double R = sqrt(pow(R0, 2) + pow(position.x, 2) + pow(position.y, 2));
                                       //    surface.height[newSurf + 2] += (R - R0);
                                       //    surface.radiusCurvature[newSurf + 2] = -R;
                                       //    surface.outerRadius[newSurf + 2] = R;
                                       //    perturbation.decenterX[newSurf + 2] = position.x + angle.x*R;
                                       //    perturbation.decenterY[newSurf + 2] = position.y + angle.y*R;
                                       //    perturbation.defocus[newSurf + 2] = position.z + angle.z*R - surface.height[newSurf + 2];
                                       //    double cosPsi(1.0), sinPsi(0.0);
                                       //    double sinPhi = -angle.x / sqrt(angle.x*angle.x + angle.y*angle.y);
                                       //    double cosPhi = angle.y / sqrt(angle.x*angle.x + angle.y*angle.y);
                                       //    double cosTheta(angle.z), sinTheta;
                                       //    if (angle.x != 0)
                                       //        sinTheta = angle.x/sinPhi;
                                       //    else if (angle.y != 0)
                                       //        sinTheta = -angle.y/cosPhi;
                                       //    else
                                       //        sinTheta = 0;
                                       //    perturbation.rotationmatrix[9*(newSurf + 2) + 0*3 + 0] = cosPsi*cosPhi - cosTheta*sinPhi*sinPsi;
                                       //    perturbation.rotationmatrix[9*(newSurf + 2) + 0*3 + 1] = cosPsi*sinPhi + cosTheta*cosPhi*sinPsi;
                                       //    perturbation.rotationmatrix[9*(newSurf + 2) + 0*3 + 2] = sinPsi*sinTheta;
                                       //    perturbation.rotationmatrix[9*(newSurf + 2) + 1*3 + 0] = -sinPsi*cosPhi - cosTheta*sinPhi*cosPsi;
                                       //    perturbation.rotationmatrix[9*(newSurf + 2) + 1*3 + 1] = -sinPsi*sinPhi + cosTheta*cosPhi*cosPsi;
                                       //    perturbation.rotationmatrix[9*(newSurf + 2) + 1*3 + 2] = cosPsi*sinTheta;
                                       //    perturbation.rotationmatrix[9*(newSurf + 2) + 2*3 + 0] = sinTheta*sinPhi;
                                       //    perturbation.rotationmatrix[9*(newSurf + 2) + 2*3 + 1] = -sinTheta*cosPhi;
                                       //    perturbation.rotationmatrix[9*(newSurf + 2) + 2*3 + 2] = cosTheta;
                                       //    for (int ii = 0; ii < 3; ii++) {
                                       //        for (int jj = 0; jj < 3; jj++) {
                                       //            perturbation.inverserotationmatrix[9*(newSurf + 2) + ii*3 + jj] =
                                       //                perturbation.rotationmatrix[9*(newSurf + 2) + jj*3 + ii];
                                       //        }
                                       //    }
                                       //    for (int surfIdx = 0; surfIdx <= nsurf + 2; surfIdx++) {
                                       //       surface.innerRadius[surfIdx] = surface.innerRadius0[surfIdx];
                                       //       surface.asphere(surfIdx, SURFACE_POINTS);
                                       //    }
                                       //}
                                       //transform(&position, &angle, newSurf + 2);
                                       //miss = findSurface(angle, position, &distance, newSurf + 2);
                                       //if (miss) {
                                       //    std::cout << "Ray "<< ray << "missed exit pupil plane\n";
                                       //}
                                       //photon.op -= distance;
                                       //transformInverse(&position, &angle, newSurf + 2);
                                       //if (eventfile) {
                                       //    Vector pt;
                                       //    vectorCopy(position, &pt);
                                       //    propagate(&pt, angle, distance);
                                       //    pEventLogging->logPhoton(pt.x, pt.y, pt.z, newSurf + 202);
                                       //}
                                       //solve line-sphere intersection analytically
                                       if (ray == 0) {
                                           //use zemax CR
                                           //position.x = 0.0;
                                           //position.y = 306.2964712316;
                                           //position.z = 4428.7341544063;
                                           cx = position.x;
                                           cy = position.y;
                                           cz = position.z;
                                           r0 = sqrt(pow(surface.outerRadius[newSurf + 2], 2) + pow(cx, 2) + pow(cy, 2));
                                           printf("chief ray position: %.9f %.9f %.9f\n",cx,cy,cz);
                                           printf("EPR: %.9f \n",surface.outerRadius[newSurf + 2]);
                                           for (int surfIdx = 0; surfIdx <= nsurf + 2; surfIdx++) {
                                              surface.innerRadius[surfIdx] = surface.innerRadius0[surfIdx];
                                              surface.asphere(surfIdx, SURFACE_POINTS);
                                           }
                                       } else if (ray == targetRay ) {
                                           printf("%g %g %g %g %g %g\n",position.x,position.y,position.z,angle.x,angle.y,angle.z);
                                       }
                                       double ocx = position.x - cx;
                                       double ocy = position.y - cy;
                                       double ocz = position.z - cz;
                                       double ocsqr = ocx*ocx + ocy*ocy + ocz*ocz;
                                       double ocproj = ocx*angle.x + ocy*angle.y + ocz*angle.z;
                                       distance = - ocproj + sqrt(ocproj*ocproj - ocsqr + r0*r0);
                                       photon.op -= distance;
                                       if (eventfile) {
                                           Vector pt;
                                           vectorCopy(position, &pt);
                                           propagate(&pt, angle, distance);
                                           pEventLogging->logPhoton(pt.x, pt.y, pt.z, newSurf + 202);
                                       }
                                       if (eventfile) pEventLogging->logPhoton(distance, 1.0, photon.op, newSurf + 402);
                                    }
                                    transform(&position, &angle, newSurf);
                                    transform(&position, &angle, newSurf + 1);
                                }
                                vectorCopy(position, &positionPrevious);
                                vectorCopy(angle, &anglePrevious);
                                detRay = 0;

                            detectorloop: while (detRay < backSourceOver) {


                                    position.x = positionPrevious.x + random_gaussian()*backSigma;
                                    position.y = positionPrevious.y + random_gaussian()*backSigma;
                                    position.z = positionPrevious.z;
                                    vectorCopy(anglePrevious, &angle);

                                    //   SILICON
                                    if (detector_on) {

                                        photon.xPos = static_cast<long>(floor(position.x*1000/pixsize + pixelsx/2));
                                        photon.yPos = static_cast<long>(floor(position.y*1000/pixsize + pixelsy/2));
                                        if (photon.xPos <= minx - activeBuffer || photon.xPos >= maxx + activeBuffer ||
                                            photon.yPos <= miny - activeBuffer || photon.yPos >= maxy + activeBuffer) {
                                            countBad(&counterLog, sourceOver, &ray);
                                            detRay++;
                                            goto detectorloop;
                                        }
                                        photon.xPosR = position.x*1000/pixsize - floor(position.x*1000/pixsize) - 0.5;
                                        photon.yPosR = position.y*1000/pixsize - floor(position.y*1000/pixsize) - 0.5;

                                        double dh;
                                        miss = getDeltaIntercept(position.x, position.y, &dh, newSurf);

                                        if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny &&
                                            photon.yPos <= maxy && contaminationmode == 1) {
                                            if (RngDouble()>(double)(*(contamination.chiptransmission +
                                                                       chip.nampx*(photon.yPos - miny) + (photon.xPos - minx)))) {
                                                countBad(&counterLog, sourceOver, &ray);
                                                detRay++;
                                                goto detectorloop;
                                            }
                                            if (*(contamination.chiplistmap + chip.nampx*(photon.yPos - miny) + (photon.xPos - minx)) != -1) {
                                                double xx = ((position.x*1000/pixsize) + pixelsx/2)*1e-3*pixsize;
                                                double yy = ((position.y*1000/pixsize) + pixelsy/2)*1e-3*pixsize;
                                                long cc = *(contamination.chiplistmap + chip.nampx*(photon.yPos - miny) + (photon.xPos - minx));
                                                if (RngDouble() > exp(-contamination.absorptionLength*contamination.chiplists[cc])) {
                                                    countBad(&counterLog, sourceOver, &ray);
                                                    detRay++;
                                                    goto detectorloop;
                                                }
                                                if (sqrt(pow(xx - contamination.chiplistx[cc], 2.0) +
                                                         pow(yy - contamination.chiplisty[cc], 2.0)) <
                                                    contamination.chiplists[cc]) {
                                                    long index;
                                                    find(contamination.henyey_greenstein, contamination.elements, RngDouble(), &index);
                                                    double mu = contamination.henyey_greenstein_mu[index];
                                                    double phi = RngDouble()*2*M_PI;
                                                    shift_mu(&angle, mu, phi);
                                                    if (mu < 0) {
                                                        countBad(&counterLog, sourceOver, &ray);
                                                        detRay++;
                                                        goto detectorloop;
                                                    }
                                                }
                                            }
                                        }


                                        miss = photonSiliconPropagate(&angle, &position, photon.wavelength, normal, dh, waveSurfaceIndex);

                                        if (miss == 0) {
                                            if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, 300);
                                        } else {
                                            if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, 302);
                                            countBad(&counterLog, sourceOver, &ray);
                                            detRay++;
                                            goto detectorloop;
                                        }

                                        miss = electronSiliconPropagate(&angle, &position);

                                        if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, 301);

                                    }


                                    photon.xPos = (long)(floor(position.x*1000/pixsize + pixelsx/2));
                                    photon.yPos = (long)(floor(position.y*1000/pixsize + pixelsy/2));

                                    if (eventfile) {
                                        if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny && photon.yPos <= maxy) {
                                            pEventLogging->logPhoton(static_cast<double>(photon.xPos),
                                                                     static_cast<double>(photon.yPos), 0.0, 303);
                                        }
                                    }

                                    if (centroidfile) {
                                        if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny && photon.yPos <= maxy) {
                                            source_xpos[source] += photon.xPos*sourceOver;
                                            source_ypos[source] += photon.yPos*sourceOver;
                                            source_photon[source] += sourceOver;
                                        }
                                    }


                                    if (opdfile) {
                                        if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny && photon.yPos <= maxy) {
                                            long xx = photon.opdx/maxr/2*(OPD_SCREEN_SIZE - 1) + OPD_SCREEN_SIZE/2.0;
                                            long yy = photon.opdy/maxr/2*(OPD_SCREEN_SIZE - 1) + OPD_SCREEN_SIZE/2.0;
                                            if (xx >= 0 && xx < OPD_SCREEN_SIZE && yy >= 0 && yy < OPD_SCREEN_SIZE) {
                                                if (ray==targetRay) {
                                                    printf("xx, yy: %.10f %i %i %.7f %.7f\n",photon.op,xx,yy,photon.opdx,photon.opdy);
                                                }
                                                *(opd + OPD_SCREEN_SIZE*xx + yy) += photon.op;
                                                *(opdcount + OPD_SCREEN_SIZE*xx + yy) += 1;
                                            }
                                        }
                                    }

                                    if (sources.type[source] < 4 && backGamma > 1.0) {
                                        Vector newPosition;
                                        vectorCopy(position, &newPosition);
                                        long long lmax;
                                        lmax = photon.sourceOver_m;
                                        photon.sourceOver_m = 1;
                                        for (long long l = 0; l < lmax; l++) {
                                            if (l > 0) {
                                                position.x = newPosition.x + random_gaussian()*backSigma/backDelta;
                                                position.y = newPosition.y + random_gaussian()*backSigma/backDelta;
                                            }
                                            photon.xPos = static_cast<long>(floor(position.x*1000/pixsize + pixelsx/2));
                                            photon.yPos = static_cast<long>(floor(position.y*1000/pixsize + pixelsy/2));
                                            if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny && photon.yPos <= maxy) {
                                                if (saturation) {
                                                    saturate(source, &largeAngle);
                                                } else {
                                                    *(chip.focal_plane + chip.nampx*(photon.yPos - miny) +
                                                      (photon.xPos - minx)) += photon.sourceOver_m;
                                                }
                                                countGood(&counterLog, photon.sourceOver_m, &ray);
                                            } else {
                                                countBad(&counterLog, photon.sourceOver_m, &ray);
                                            }
                                        }
                                        photon.sourceOver_m = lmax;
                                    } else {
                                        if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny && photon.yPos <= maxy) {
                                            if (saturation) {
                                                saturate(source, &largeAngle);
                                            } else {
                                                *(chip.focal_plane + chip.nampx*(photon.yPos - miny) + (photon.xPos - minx)) += photon.sourceOver_m;
                                            }
                                        } else {
                                            countBad(&counterLog, sourceOver, &ray);
                                            detRay++;
                                            goto detectorloop;
                                        }
                                    }


                                    if (throughputfile) addThroughput(&throughputLog, nsurf, waveIndex, photon.sourceOver_m);
                                    detRay++;
                                    if (sources.type[source] >= 4 || (sources.type[source] < 4 && backGamma <= 1.0)) {
                                        countGood(&counterLog, photon.sourceOver_m, &ray);
                                    }
                                }
                                break;

                            }
                            //if ( ray == 0 ) printf("D: surf %d\n",newSurf);
                            //if ( ray == 0 ) printf("%g %g %g %g %g %g\n",position.x,position.y,position.z,angle.x,angle.y,angle.z);
                        }
                    }
                    sourceCounter++;
                }
            }
        }

        if (sourceType >= 0  && sourceType < 50) sprintf(tempstring, "Dome Light         ");
        if (sourceType >= 50 && sourceType < 100) sprintf(tempstring, "Dark Sky           ");
        if (sourceType >= 100 && sourceType < 150) sprintf(tempstring, "Moon               ");
        if (sourceType >= 150 && sourceType < 200) sprintf(tempstring, "Zodiacal Light     ");
        if (sourceType == 200) sprintf(tempstring, "Astrophysical m>40 ");
        if (sourceType >= 201 && sourceType < 242) sprintf(tempstring, "Astrophysical m=%2d ", 241 - sourceType);
        if (sourceType == 242) sprintf(tempstring, "Astrophysical m<0  ");
        if (sourceCounter > 0) counterCheck(&counterLog, sourceCounter, tempstring);
    }

    // COSMIC RAYS
    if (checkpointcount == checkpointtotal) {
        detRay = 0;
        ray = 0;
        cosmicRays(&detRay);
        if (detRay > 0) {
            for (long i = 0; i < detRay; i++) {
                countGood(&counterLog, 1, &ray);
            }
            sprintf(tempstring, "Cosmic Rays        ");
            counterCheck(&counterLog, detRay, tempstring);
        }
    }

    // OUTPUT DATA
    eventfile=writeEvent;
    if (checkpointcount == checkpointtotal) writeImageFile();
    if (opdfile) writeOPD();
    if (checkpointcount !=  checkpointtotal) writeCheckpoint(checkpointcount);
    if (centroidfile) writeCentroidFile(outputdir, outputfilename, source_photon, source_xpos, source_ypos, sources.id, nsource);
    if (throughputfile) writeThroughputFile(outputdir, outputfilename, &throughputLog, nsurf);
    if (eventfile) pEventLogging->eventFileClose();
    return(0);

}
