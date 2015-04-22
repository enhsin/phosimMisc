///
/// @package phosim
/// @file photonmanipulate.cpp
/// @brief photon manipulation routines (part of Image class)
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

int Image::domeSeeing(Vector *angle) {

    double phi, r;

    phi = 2*M_PI*RngDouble();
    r = sqrt(domeseeing*domeseeing + toypsf*toypsf)*ARCSEC/2.35482*random_gaussian();

    angle->x = angle->x + r*cos(phi);
    angle->y = angle->y + r*sin(phi);
    angle->z = smallAnglePupilNormalize(angle->x, angle->y);

    return(0);

}

int Image::tracking(Vector *angle, double time) {

    double rindex, vxp, vyp;
    long index;

    index = find_linear(perturbation.jittertime, trackinglines, time, &rindex);
    vxp = angle->x - (cos(spiderangle)*perturbation.jitterele[index] +
                      sin(spiderangle)*perturbation.jitterazi[index])*ARCSEC;
    vyp = angle->y + ((-1.0)*(sin(spiderangle)*perturbation.jitterele[index]) +
                      cos(spiderangle)*perturbation.jitterazi[index])*ARCSEC;
    angle->y = vyp*cos(perturbation.jitterrot[index]*ARCSEC) + (-1.0)*vxp*sin(perturbation.jitterrot[index]*ARCSEC);
    angle->x = vyp*sin(perturbation.jitterrot[index]*ARCSEC) + (vxp)*cos(perturbation.jitterrot[index]*ARCSEC);
    angle->z = smallAnglePupilNormalize(angle->x, angle->y);

    return(0);
}

int Image::largeAngleScattering(Vector *largeAngle) {

    long index;
    double r, phi;

    r = RngDouble();
    if (r < miescatter_scat) {
        r = RngDouble();
        find(perturbation.miescatter, 10000, r, &index);
        r = ((double)(index))/10000.0*0.1*DEGREE;
        phi = 2.0*M_PI*RngDouble();
        largeAngle->x += r*cos(phi);
        largeAngle->y += r*sin(phi);
    }

    return(0);

}

int Image::secondKick(Vector *largeAngle) {

    long index;
    double r, phi;

    r = RngDouble();
    find(screen.hffunc, 10000, r, &index);
    r = ((double)(index))*0.5*1e-3/(SCREEN_SIZE*screen.fine_sizeperpixel)*photon.wavelengthFactor;
    phi = 2*M_PI*RngDouble();
    largeAngle->x += r*cos(phi);
    largeAngle->y += r*sin(phi);

    return(0);

}


int Image::diffraction(Vector *position, Vector angle, Vector *largeAngle) {

    double distance;
    double xl, yl, zorig, r;
    int i;
    double diffdist, adiffdist;
    int signv = 0;
    double diffx = 0.0, diffy = 0.0, mindiffdist;

    mindiffdist = 1e30;
    zorig = position->z;
    double cosShiftedAngle = cos(photon.shiftedAngle);
    double sinShiftedAngle = sin(photon.shiftedAngle);

    for (i = 0; i < obstruction.nspid; i++) {

        distance = (obstruction.spider_height[i] - (position->z))/angle.z;
        propagate(position, angle, distance);
        xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
        yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;

        if (obstruction.spider_type[i] == 1) {

            if (obstruction.spider_angle[i] != 0) {
                distance = (obstruction.spider_height[i] +
                            (fabs(yl) - obstruction.spider_reference[i])*sin(DEGREE*obstruction.spider_angle[i]) -
                            (position->z))/angle.z;
                propagate(position, angle, distance);
                xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
                yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;
            }

            diffdist = xl - obstruction.spider_center[i];
            if (diffdist > 0) {
                signv = 1;
            } else {
                signv = -1;
            }
            adiffdist = fabs(diffdist);
            if (adiffdist < obstruction.spider_width[i]) return(1);
            if (adiffdist - obstruction.spider_width[i] < (mindiffdist)) {
                mindiffdist = fabs(diffdist) - obstruction.spider_width[i];
                diffx = cosShiftedAngle*diffdist/adiffdist;
                diffy = -sinShiftedAngle*diffdist/adiffdist;
            }

            distance = (obstruction.spider_height[i] - obstruction.spider_depth[i] +
                        (fabs(yl) - obstruction.spider_reference[i])*sin(DEGREE*obstruction.spider_angle[i]) -
                        (position->z))/angle.z;
            propagate(position, angle, distance);
            xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
            yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;

            diffdist = xl - obstruction.spider_center[i];
            if (diffdist > 0 && signv == -1) return(1);
            if (diffdist < 0 && signv == 1) return(1);
            adiffdist = fabs(diffdist);
            if (adiffdist < obstruction.spider_width[i]) return(1);
            if (adiffdist - obstruction.spider_width[i] < (mindiffdist)) {
                mindiffdist = fabs(diffdist) - obstruction.spider_width[i];
                diffx = cosShiftedAngle*diffdist/adiffdist;
                diffy = -sinShiftedAngle*diffdist/adiffdist;
            }

        }

        if (obstruction.spider_type[i] == 2) {

            if (obstruction.spider_angle[i] != 0) {
                distance = (obstruction.spider_height[i] +
                            (fabs(xl)-obstruction.spider_reference[i])*sin(DEGREE*obstruction.spider_angle[i]) -
                            (position->z))/angle.z;
                propagate(position, angle, distance);
                xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
                yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;
            }

            diffdist = yl - obstruction.spider_center[i];
            if (diffdist > 0) {
                signv = 1;
            } else {
                signv = -1;
            }
            adiffdist = fabs(diffdist);
            if (adiffdist < obstruction.spider_width[i]) return(1);
            if (adiffdist - obstruction.spider_width[i] < (mindiffdist)) {
                mindiffdist = fabs(diffdist) - obstruction.spider_width[i];
                diffx = sinShiftedAngle*diffdist/adiffdist;
                diffy = cosShiftedAngle*diffdist/adiffdist;
            }

            distance = (obstruction.spider_height[i] - obstruction.spider_depth[i] +
                        (fabs(xl)-obstruction.spider_reference[i])*sin(DEGREE*obstruction.spider_angle[i]) -
                        (position->z))/angle.z;
            propagate(position, angle, distance);
            xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
            yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;

            diffdist = yl - obstruction.spider_center[i];
            if (diffdist > 0 && signv == -1) return(1);
            if (diffdist < 0 && signv == 1) return(1);
            adiffdist = fabs(diffdist);
            if (adiffdist < obstruction.spider_width[i]) return(1);
            if (adiffdist - obstruction.spider_width[i] < (mindiffdist)) {
                mindiffdist = fabs(diffdist) - obstruction.spider_width[i];
                diffx = sinShiftedAngle*diffdist/adiffdist;
                diffy = cosShiftedAngle*diffdist/adiffdist;
            }
        }
    }

    distance = (zorig - (position->z))/angle.z;
    propagate(position, angle, distance);

    if (obstruction.pupil == 1) {
        r = sqrt((position->x)*(position->x) + (position->y)*(position->y));
        diffdist = surface.outerRadius[0] - r;
        if (diffdist < 0) return(1);
        if (diffdist < (mindiffdist) ) {
            (mindiffdist) = diffdist;
            (diffx) = (position->x)/r;
            (diffy) = (position->y)/r;
        }
        diffdist = r - surface.innerRadius[0];
        if (diffdist < 0) return(1);
        if (diffdist < (mindiffdist) ) {
            (mindiffdist) = diffdist;
            (diffx) = (position->x)/r;
            (diffy) = (position->y)/r;
        }
    }

    r = (photon.wavelength/1000.0)/(4*M_PI*mindiffdist);
    if (diffraction_on != 5) {
        largeAngle->x += r*diffx;
        largeAngle->y += r*diffy;
    }

    return(0);

}

double Image::airIndexRefraction() {

    double airRefraction = 64.328 + 29498.1/(146 - 1/photon.wavelength/photon.wavelength) + 255.4/(41 - 1/photon.wavelength/photon.wavelength);
    airRefraction = airRefraction*pressure*(1 + (1.049 - 0.0157*temperature)*1e-6*pressure)/720.883/(1 + 0.003661*temperature);
    airRefraction = airRefraction - ((0.0624 - 0.000680/photon.wavelength/photon.wavelength)/(1 + 0.003661*temperature)*water_pressure);
    if (airrefraction) {
        return(airRefraction);
    } else {
        return(0.0);
    }

}

int Image::atmosphericDispersion (Vector *angle) {

    double dx, dy, adcx = 0.0, adcy = 0.0;

    if (atmosphericdispcenter) {
        dx = zenith*sin(photon.shiftedAngle);
        dy = zenith*cos(photon.shiftedAngle);
        adcx = tan(sqrt(dx*dx + dy*dy))*dx/sqrt(dx*dx + dy*dy)*air.air_refraction_adc/1e6;
        adcy = tan(sqrt(dx*dx + dy*dy))*dy/sqrt(dx*dx + dy*dy)*air.air_refraction_adc/1e6;
        if (zenith==0.0) {
            adcx = 0.0;
            adcy = 0.0;
        }
    }
    if (atmospheric_dispersion) {
        dx = -angle->x + zenith*sin(photon.shiftedAngle);
        dy = -angle->y + zenith*cos(photon.shiftedAngle);
        angle->x = angle->x + tan(sqrt(dx*dx + dy*dy))*dx/sqrt(dx*dx + dy*dy)*photon.airRefraction/1e6;
        angle->y = angle->y + tan(sqrt(dx*dx + dy*dy))*dy/sqrt(dx*dx + dy*dy)*photon.airRefraction/1e6;
        if (atmosphericdispcenter) {
            angle->x = angle->x - adcx;
            angle->y = angle->y - adcy;
        }
        angle->z = smallAnglePupilNormalize(angle->x, angle->y);
    }
    return(0);

}


int Image::samplePupil (Vector *position, long long ray) {

    double r, phi, rindex;
    long index;

    if (aperturemode == 1) {
        double x, y;
        x = (((double)((ray % 40000) % 200))/200.0)*2.0*maxr - maxr;
        y = (((double)((ray % 40000) / 200))/200.0)*2.0*maxr - maxr;
        r = sqrt(x*x + y*y);
        phi = atan2(y, x);
        if (r < minr || r > maxr) return(1);
    } else if (aperturemode == 2) {
        if ( ray == 0 ) {
            r = 1e-14;
            phi = 0.0;
        } else {
            double x, y;
            long OPDSIZE2 = OPD_SCREEN_SIZE*OPD_SCREEN_SIZE;
            x = (((double)((ray % OPDSIZE2 ) % OPD_SCREEN_SIZE))/(OPD_SCREEN_SIZE - 1))*2.0*maxr - maxr;
            y = (((double)((ray % OPDSIZE2 ) / OPD_SCREEN_SIZE))/(OPD_SCREEN_SIZE - 1))*2.0*maxr - maxr;
            r = sqrt(x*x + y*y);
            phi = atan2(y, x);
            if (r < minr || r > maxr) return(1);
        }
    } else {
        r = sqrt(RngDouble()*(maxr*maxr - minr*minr) + minr*minr);
        phi = RngDouble()*2*M_PI;
    }

    position->x = r*cos(phi);
    position->y = r*sin(phi);
    if (finiteDistance == 0.0 && aperturemode != 2) {
        index = find_linear(&surface.radius[0], SURFACE_POINTS, r, &rindex);
        position->z = interpolate_linear(&surface.profile[0], index, rindex);
    } else {
        position->z = surface.height[0];
    }

    return(0);

}


int Image::transmissionCheck(double transmission, long surfaceIndex, long waveSurfaceIndex) {

    double randNum;

    (photon.counter)++;
    if (transmission > dynamicTransmission[waveSurfaceIndex+surfaceIndex]) {
        dynamicTransmission[waveSurfaceIndex+surfaceIndex] = transmission;
    }
    if (photon.counter <= photon.maxcounter) {
        randNum = saveRand[photon.counter];
    } else {
        randNum = RngDouble();
    }
    if (randNum > transmission) {
        return(1);
    } else {
        return(0);
    }

}

int Image::transmissionPreCheck(long surfaceIndex, long waveSurfaceIndex) {

    (photon.counter)++;
    saveRand[photon.counter] = RngDouble();
    if (saveRand[photon.counter] > (fabs(dynamicTransmission[waveSurfaceIndex+surfaceIndex]) + transtol)) {
        return(1);
    } else {
        return(0);
    }

}

double Image::surfaceCoating (double wavelength, Vector angle, Vector normal, long newSurf, double *reflection) {

    double filterAngle, rindex, crindex;
    long index, cindex;

    if (surface.surfacecoating[newSurf] != 0 && coatingmode == 1) {
        if (coating.angleNumber[newSurf] > 1) {
            Vector tempangle;
            vectorCopy(angle, &tempangle);
            refract(&tempangle, normal, photon.ncurr, 1.0);
            double arg = fabs(normal.x*tempangle.x + normal.y*tempangle.y + normal.z*tempangle.z);
            if (arg > 1) arg = 1.0;
            filterAngle = acos(arg)/DEGREE;
            index = find_linear(coating.wavelength[newSurf], coating.wavelengthNumber[newSurf], wavelength, &rindex);
            cindex = find_linear(coating.angle[newSurf], coating.angleNumber[newSurf], filterAngle, &crindex);
            if (crindex - cindex < 0.0) crindex = static_cast<double>(cindex);
            if (crindex - cindex > 1.0) crindex = static_cast<double>(cindex + 1);
            if (rindex - index < 0.0) rindex = static_cast<double>(index);
            if (rindex - index > 1.0) rindex = static_cast<double>(index + 1);
            *reflection = interpolate_bilinear(coating.reflection[newSurf], coating.wavelengthNumber[newSurf], cindex, crindex, index, rindex);
            return(interpolate_bilinear(coating.transmission[newSurf], coating.wavelengthNumber[newSurf], cindex, crindex, index, rindex));
        } else {
            index = find_linear(coating.wavelength[newSurf], coating.wavelengthNumber[newSurf], wavelength, &rindex);
            if (rindex - index < 0.0) rindex = static_cast<double>(index);
            if (rindex - index > 1.0) rindex = static_cast<double>(index + 1);
            *reflection = interpolate_linear(coating.reflection[newSurf], index, rindex);
            return(interpolate_linear(coating.transmission[newSurf], index, rindex));
        }
    } else {
        return(1.0);
    }

}

void Image::newRefractionIndex(long surfaceIndex) {

    long index, newMedium;
    double rindex;

    photon.nprev = photon.ncurr;
    if (photon.direction == 1) {
        newMedium = surfaceIndex;
    } else {
        newMedium = surfaceIndex - 1;
    }
    if (newMedium >= 0) {
        if (surface.surfacemed[newMedium] == 0) {
            photon.ncurr = 1.0;
        } else if (surface.surfacemed[newMedium] == 2) {
            photon.ncurr = 1.0 + photon.airRefraction/1e6;
        } else {
            index = find_linear(medium.indexRefractionWavelength[newMedium], medium.indexRefractionNumber[newMedium], photon.wavelength, &rindex);
            photon.ncurr = interpolate_linear(medium.indexRefraction[newMedium], index, rindex);
            if (opdfile) photon.ncurr = 1.4538548621; //to compare with zemax
        }
    } else {
        photon.ncurr = 1.0 + photon.airRefraction/1e6;
    }

}


void Image::atmospherePropagate(Vector *position, Vector angle, long layer, int mode) {

    if (layer == -1) {

        propagate(position, angle, (1e6*70.0 - position->z)/(angle.z));
        photon.xporig = position->x;
        photon.yporig = position->y;
        photon.zporig = position->z;

        if (mode == 2) {
            if (fabs(photon.time - photon.prtime) > screentol) {
                for (long k = 0; k < SCREEN_SIZE; k++) {
                    for (long l = 0; l < SCREEN_SIZE; l++) {
                        *(screen.phasescreen + k*SCREEN_SIZE + l) = 0;
                    }
                }
            }
        }

    } else {

        double distance = (1e6*height[layer] - position->z)/(angle.z);
        propagate(position, angle, distance);


    }

}

void Image::atmosphereIntercept(Vector *position, long layer) {

    double rindex;
    long index;
    double wx = wind[layer]*1.0e3*(photon.absoluteTime)*cos(winddir[layer]*DEGREE - azimuth);
    double wy = wind[layer]*1.0e3*(photon.absoluteTime)*sin(winddir[layer]*DEGREE - azimuth);

    // wind blur
    index = find_linear(perturbation.jittertime, trackinglines, photon.absoluteTime, &rindex);

    photon.windx = wx - wy*screen.jitterwind[index]*DEGREE;
    photon.windy = wy + wx*screen.jitterwind[index]*DEGREE;

    photon.lindex = SCREEN_SIZE*SCREEN_SIZE*layer;

    photon.xpos = position->x + xtelloc + photon.windx;
    photon.ypos = position->y + ytelloc + photon.windy;

    find_linear_wrap(photon.xpos, screen.large_sizeperpixel, SCREEN_SIZE, &(photon.indexlx0), &(photon.indexlx1), &(photon.dlx));
    find_linear_wrap(photon.ypos, screen.large_sizeperpixel, SCREEN_SIZE, &(photon.indexly0), &(photon.indexly1), &(photon.dly));
    find_linear_wrap(photon.xpos, screen.coarse_sizeperpixel, SCREEN_SIZE, &(photon.indexcx0), &(photon.indexcx1), &(photon.dcx));
    find_linear_wrap(photon.ypos, screen.coarse_sizeperpixel, SCREEN_SIZE, &(photon.indexcy0), &(photon.indexcy1), &(photon.dcy));
    find_linear_wrap(photon.xpos, screen.medium_sizeperpixel, SCREEN_SIZE, &(photon.indexmx0), &(photon.indexmx1), &(photon.dmx));
    find_linear_wrap(photon.ypos, screen.medium_sizeperpixel, SCREEN_SIZE, &(photon.indexmy0), &(photon.indexmy1), &(photon.dmy));
    find_linear_wrap(photon.xpos, screen.fine_sizeperpixel, SCREEN_SIZE, &(photon.indexfx0), &(photon.indexfx1), &(photon.dfx));
    find_linear_wrap(photon.ypos, screen.fine_sizeperpixel, SCREEN_SIZE, &(photon.indexfy0), &(photon.indexfy1), &(photon.dfy));

}

void Image::atmosphereRefraction (Vector *angle, long layer, int mode) {

    double scaleOuter;

    if (mode == 1 || mode == 5) {
        scaleOuter = photon.wavelengthFactor*ARCSEC*screen.secondKickSize;
    } else {
        scaleOuter = photon.wavelengthFactor*ARCSEC;
    }

    if (mode <= 1 || mode >= 4) {

        if (atmdebug == 0) {

            (angle->x) += (interpolate_bilinear_float_wrap(screen.seex_large + photon.lindex, SCREEN_SIZE, photon.indexlx0, photon.indexlx1,
                                                         photon.dlx, photon.indexly0, photon.indexly1, photon.dly) +
                         interpolate_bilinear_float_wrap(screen.seex_coarse + photon.lindex, SCREEN_SIZE, photon.indexcx0, photon.indexcx1,
                                                         photon.dcx, photon.indexcy0, photon.indexcy1, photon.dcy) +
                         interpolate_bilinear_float_wrap(screen.seex_medium + photon.lindex, SCREEN_SIZE, photon.indexmx0, photon.indexmx1,
                                                         photon.dmx, photon.indexmy0, photon.indexmy1, photon.dmy))*scaleOuter;

            (angle->y) += (interpolate_bilinear_float_wrap(screen.seey_large + photon.lindex, SCREEN_SIZE, photon.indexlx0, photon.indexlx1,
                                                         photon.dlx, photon.indexly0, photon.indexly1, photon.dly) +
                         interpolate_bilinear_float_wrap(screen.seey_coarse + photon.lindex, SCREEN_SIZE, photon.indexcx0, photon.indexcx1,
                                                         photon.dcx, photon.indexcy0, photon.indexcy1, photon.dcy) +
                         interpolate_bilinear_float_wrap(screen.seey_medium + photon.lindex, SCREEN_SIZE, photon.indexmx0, photon.indexmx1,
                                                         photon.dmx, photon.indexmy0, photon.indexmy1, photon.dmy))*scaleOuter;

            angle->z = smallAnglePupilNormalize(angle->x, angle->y);

        } else {

            photon.indexlx0 = static_cast<long>(floor((photon.indexlx0)/large_grid)*large_grid);
            photon.indexly0 = static_cast<long>(floor((photon.indexly0)/large_grid)*large_grid);
            photon.indexcx0 = static_cast<long>(floor((photon.indexcx0)/coarse_grid)*coarse_grid);
            photon.indexcy0 = static_cast<long>(floor((photon.indexcy0)/coarse_grid)*coarse_grid);
            photon.indexmx0 = static_cast<long>(floor((photon.indexmx0)/medium_grid)*medium_grid);
            photon.indexmy0 = static_cast<long>(floor((photon.indexmy0)/medium_grid)*medium_grid);
            photon.indexfx0 = static_cast<long>(floor((photon.indexfx0)/fine_grid)*fine_grid);
            photon.indexfy0 = static_cast<long>(floor((photon.indexfy0)/fine_grid)*fine_grid);
            photon.indexlx1 = static_cast<long>(floor((photon.indexlx1)/large_grid)*large_grid);
            photon.indexly1 = static_cast<long>(floor((photon.indexly1)/large_grid)*large_grid);
            photon.indexcx1 = static_cast<long>(floor((photon.indexcx1)/coarse_grid)*coarse_grid);
            photon.indexcy1 = static_cast<long>(floor((photon.indexcy1)/coarse_grid)*coarse_grid);
            photon.indexmx1 = static_cast<long>(floor((photon.indexmx1)/medium_grid)*medium_grid);
            photon.indexmy1 = static_cast<long>(floor((photon.indexmy1)/medium_grid)*medium_grid);
            photon.indexfx1 = static_cast<long>(floor((photon.indexfx1)/fine_grid)*fine_grid);
            photon.indexfy1 = static_cast<long>(floor((photon.indexfy1)/fine_grid)*fine_grid);

            (angle->x) += (interpolate_bilinear_float_wrap(screen.seex_large + photon.lindex, SCREEN_SIZE, photon.indexlx0, photon.indexlx1,
                                                         photon.dlx, photon.indexly0, photon.indexly1, photon.dly)*large_scale +
                         interpolate_bilinear_float_wrap(screen.seex_coarse + photon.lindex, SCREEN_SIZE, photon.indexcx0, photon.indexcx1,
                                                         photon.dcx, photon.indexcy0, photon.indexcy1, photon.dcy)*coarse_scale +
                         interpolate_bilinear_float_wrap(screen.seex_medium + photon.lindex, SCREEN_SIZE, photon.indexmx0, photon.indexmx1,
                                                         photon.dmx, photon.indexmy0, photon.indexmy1, photon.dmy)*medium_scale)*scaleOuter;

            (angle->y) += (interpolate_bilinear_float_wrap(screen.seey_large + photon.lindex, SCREEN_SIZE, photon.indexlx0, photon.indexlx1,
                                                         photon.dlx, photon.indexly0, photon.indexly1, photon.dly)*large_scale +
                         interpolate_bilinear_float_wrap(screen.seey_coarse + photon.lindex, SCREEN_SIZE, photon.indexcx0, photon.indexcx1,
                                                         photon.dcx, photon.indexcy0, photon.indexcy1, photon.dcy)*coarse_scale +
                         interpolate_bilinear_float_wrap(screen.seey_medium + photon.lindex, SCREEN_SIZE, photon.indexmx0, photon.indexmx1,
                                                         photon.dmx, photon.indexmy0, photon.indexmy1, photon.dmy)*medium_scale)*scaleOuter;

            angle->z = smallAnglePupilNormalize(angle->x, angle->y);


        }

    } else {
        if (fabs(photon.time - photon.prtime) > screentol) {
            double randomi = RngDouble();
            double randomj = RngDouble();
            for (long i = 0; i < SCREEN_SIZE; i++) {
                for (long j = 0; j < SCREEN_SIZE; j++) {

                    find_linear_wrap(photon.xpos - photon.xp + (i + randomi - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.large_sizeperpixel, SCREEN_SIZE, &(photon.indexlx0), &(photon.indexlx1), &(photon.dlx));
                    find_linear_wrap(photon.xpos - photon.xp + (i + randomi - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.coarse_sizeperpixel, SCREEN_SIZE, &(photon.indexcx0), &(photon.indexcx1), &(photon.dcx));
                    find_linear_wrap(photon.xpos - photon.xp + (i + randomi - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.medium_sizeperpixel, SCREEN_SIZE, &(photon.indexmx0), &(photon.indexmx1), &(photon.dmx));
                    find_linear_wrap(photon.xpos - photon.xp + (i + randomi - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.fine_sizeperpixel,  SCREEN_SIZE, &(photon.indexfx0), &(photon.indexfx1), &(photon.dfx));
                    find_linear_wrap(photon.ypos - photon.yp + (j + randomj - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.large_sizeperpixel, SCREEN_SIZE, &(photon.indexly0), &(photon.indexly1), &(photon.dly));
                    find_linear_wrap(photon.ypos - photon.yp + (j + randomj - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.coarse_sizeperpixel, SCREEN_SIZE, &(photon.indexcy0), &(photon.indexcy1), &(photon.dcy));
                    find_linear_wrap(photon.ypos - photon.yp + (j + randomj - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.medium_sizeperpixel, SCREEN_SIZE, &(photon.indexmy0), &(photon.indexmy1), &(photon.dmy));
                    find_linear_wrap(photon.ypos - photon.yp + (j + randomj - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.fine_sizeperpixel,  SCREEN_SIZE,  &(photon.indexfy0),  &(photon.indexfy1), &(photon.dfy));

                    if (atmdebug == 0) {
                        if (mode == 2) {

                            *(screen.phasescreen + i*SCREEN_SIZE + j) +=
                                (interpolate_bilinear_float_wrap(screen.phase_large + photon.lindex, SCREEN_SIZE, photon.indexlx0, photon.indexlx1,
                                                                 photon.dlx, photon.indexly0, photon.indexly1, photon.dly)+
                                 interpolate_bilinear_float_wrap(screen.phase_coarse + photon.lindex,SCREEN_SIZE, photon.indexcx0, photon.indexcx1,
                                                                 photon.dcx, photon.indexcy0, photon.indexcy1, photon.dcy)+
                                 interpolate_bilinear_float_wrap(screen.phase_medium + photon.lindex,SCREEN_SIZE, photon.indexmx0, photon.indexmx1,
                                                                 photon.dmx, photon.indexmy0, photon.indexmy1, photon.dmy)+
                                 interpolate_bilinear_float_wrap(screen.phase_fine + photon.lindex  ,SCREEN_SIZE, photon.indexfx0, photon.indexfx1,
                                                                 photon.dfx, photon.indexfy0, photon.indexfy1, photon.dfy))*
                                screen.phase_norm[layer]*(seefactor[layer]/(totalseeing*pow(1/cos(zenith), 0.6)/2.35));

                        } else {

                            *(screen.phasescreen + i*SCREEN_SIZE + j) +=
                                (interpolate_bilinear_float_wrap(screen.phaseh_medium + photon.lindex, SCREEN_SIZE, photon.indexmx0, photon.indexmx1,
                                                                 photon.dmx, photon.indexmy0, photon.indexmy1, photon.dmy)+
                                 interpolate_bilinear_float_wrap(screen.phaseh_fine + photon.lindex, SCREEN_SIZE, photon.indexfx0, photon.indexfx1,
                                                                 photon.dfx,  photon.indexfy0, photon.indexfy1, photon.dfy))*
                                screen.phase_norm[layer]*(seefactor[layer]/(totalseeing*pow(1/cos(zenith), 0.6)/2.35));

                        }

                    } else {

                        photon.indexlx0 = static_cast<long>(floor((photon.indexlx0)/large_grid)*large_grid);
                        photon.indexly0 = static_cast<long>(floor((photon.indexly0)/large_grid)*large_grid);
                        photon.indexcx0 = static_cast<long>(floor((photon.indexcx0)/coarse_grid)*coarse_grid);
                        photon.indexcy0 = static_cast<long>(floor((photon.indexcy0)/coarse_grid)*coarse_grid);
                        photon.indexmx0 = static_cast<long>(floor((photon.indexmx0)/medium_grid)*medium_grid);
                        photon.indexmy0 = static_cast<long>(floor((photon.indexmy0)/medium_grid)*medium_grid);
                        photon.indexfx0 = static_cast<long>(floor((photon.indexfx0)/fine_grid)*fine_grid);
                        photon.indexfy0 = static_cast<long>(floor((photon.indexfy0)/fine_grid)*fine_grid);
                        photon.indexlx1 = static_cast<long>(floor((photon.indexlx1)/large_grid)*large_grid);
                        photon.indexly1 = static_cast<long>(floor((photon.indexly1)/large_grid)*large_grid);
                        photon.indexcx1 = static_cast<long>(floor((photon.indexcx1)/coarse_grid)*coarse_grid);
                        photon.indexcy1 = static_cast<long>(floor((photon.indexcy1)/coarse_grid)*coarse_grid);
                        photon.indexmx1 = static_cast<long>(floor((photon.indexmx1)/medium_grid)*medium_grid);
                        photon.indexmy1 = static_cast<long>(floor((photon.indexmy1)/medium_grid)*medium_grid);
                        photon.indexfx1 = static_cast<long>(floor((photon.indexfx1)/fine_grid)*fine_grid);
                        photon.indexfy1 = static_cast<long>(floor((photon.indexfy1)/fine_grid)*fine_grid);

                        if (mode == 2) {

                            *(screen.phasescreen + i*SCREEN_SIZE + j) +=
                                (interpolate_bilinear_float_wrap(screen.phase_large + photon.lindex, SCREEN_SIZE, photon.indexlx0, photon.indexlx1,
                                                                 photon.dlx, photon.indexly0, photon.indexly1, photon.dly)*large_scale +
                                 interpolate_bilinear_float_wrap(screen.phase_coarse + photon.lindex, SCREEN_SIZE, photon.indexcx0, photon.indexcx1,
                                                                 photon.dcx, photon.indexcy0, photon.indexcy1, photon.dcy)*coarse_scale +
                                 interpolate_bilinear_float_wrap(screen.phase_medium + photon.lindex, SCREEN_SIZE, photon.indexmx0, photon.indexmx1,
                                                                 photon.dmx, photon.indexmy0, photon.indexmy1, photon.dmy)*medium_scale +
                                 interpolate_bilinear_float_wrap(screen.phase_fine + photon.lindex, SCREEN_SIZE, photon.indexfx0, photon.indexfx1,
                                                                 photon.dfx, photon.indexfy0, photon.indexfy1, photon.dfy)*fine_scale)*
                                screen.phase_norm[layer]*(seefactor[layer]/(totalseeing*pow(1/cos(zenith), 0.6)/2.35));

                        } else {
                            *(screen.phasescreen + i*SCREEN_SIZE + j) +=
                                (interpolate_bilinear_float_wrap(screen.phaseh_medium + photon.lindex, SCREEN_SIZE, photon.indexmx0, photon.indexmx1,
                                                                 photon.dmx, photon.indexmy0, photon.indexmy1, photon.dmy)*medium_scale +
                                 interpolate_bilinear_float_wrap(screen.phaseh_fine + photon.lindex, SCREEN_SIZE, photon.indexfx0, photon.indexfx1,
                                                                 photon.dfx, photon.indexfy0, photon.indexfy1, photon.dfy)*fine_scale)*
                                screen.phase_norm[layer]*(seefactor[layer]/(totalseeing*pow(1/cos(zenith), 0.6)/2.35));

                        }


                    }
                }
            }
        }

    }

}


double Image::cloudOpacity (long layer) {

    double transmission;

    if (cloudmean[layer] != 0 || cloudvary[layer] != 0) {
        transmission = pow(10.0, -0.4*(cloudmean[layer] +
                                       cloudvary[layer]*((double)(*(screen.cloud[layer] + photon.indexcx0*SCREEN_SIZE + photon.indexcy0)))));
        if (transmission > 1.0) transmission = 1.0;
        return(transmission);
    } else {
        return(1.0);
    }

}


double Image::atmosphereOpacity (Vector angle, long layer) {

    double densityScale, dvx, dvy, airmassl;
    double rindex;

    if (layer == -1) {

        dvx = -angle.x + zenith*sin(photon.shiftedAngle);
        dvy = -angle.y + zenith*cos(photon.shiftedAngle);
        photon.dvr = sqrt(dvx*dvx + dvy*dvy);
        if (photon.dvr == 0.0) {
            airmassl = 1.0;
        } else {
            airmassl = RADIUS_EARTH*(sqrt(pow(1 + 70.0/RADIUS_EARTH, 2.0)/sin(photon.dvr)/sin(photon.dvr) - 1) -
                                     sqrt(pow(1 + height[layer + 1]/RADIUS_EARTH, 2.0)/sin(photon.dvr)/sin(photon.dvr) - 1))*
                sin(photon.dvr)/(70.0 - height[layer + 1]);
        }

        photon.oindex = find_linear(air.tauWavelength, 180001, photon.wavelength, &rindex);
        return(exp(-(*(air.tau[0] + photon.oindex))*airmassl/air.airmassLayer[layer + 1]));


    } else {

        if (photon.dvr == 0.0) {
            airmassl = 1.0;
        } else {
            airmassl = RADIUS_EARTH*(sqrt(pow(1 + height[layer]/RADIUS_EARTH, 2.0)/sin(photon.dvr)/sin(photon.dvr) - 1) -
                                     sqrt(pow(1 + height[layer + 1]/RADIUS_EARTH, 2.0)/sin(photon.dvr)/sin(photon.dvr) - 1))*
                sin(photon.dvr)/(height[layer] - height[layer + 1]);
        }

        if (densityfluctuation[layer] != 0.0) {
            densityScale = ((double)(*(screen.phase_large + photon.lindex + photon.indexlx0*SCREEN_SIZE + photon.indexly0)+
                                     *(screen.phase_coarse + photon.lindex + photon.indexcx0*SCREEN_SIZE + photon.indexcy0)+
                                     *(screen.phase_medium + photon.lindex + photon.indexmx0*SCREEN_SIZE + photon.indexmy0)+
                                     *(screen.phase_fine + photon.lindex + photon.indexfx0*SCREEN_SIZE + photon.indexfy0)))*densityfluctuation[layer]+densitymean[layer];
        } else {
            densityScale = densitymean[layer];
        }

        return(exp(-(*(air.tau[layer + 1] + photon.oindex))*densityScale*airmassl/air.airmassLayer[layer + 1]));

    }

}


// transform to optic frame
void Image::transform( Vector *position, Vector *angle, long surfaceIndex) {

    double vprime1, vprime2, vprime3;

    // defocus
    position->z  =  position->z - perturbation.defocus[surfaceIndex] - surface.height[surfaceIndex];

    // decenter
    position->x = position->x - perturbation.decenterX[surfaceIndex] - surface.centerx[surfaceIndex];
    position->y = position->y - perturbation.decenterY[surfaceIndex] - surface.centery[surfaceIndex];

    // angles
    vprime1 = (*(perturbation.rotationmatrix + 9*surfaceIndex+0*3 + 0))*(angle->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 1))*(angle->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 2))*(angle->z);
    vprime2 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 0))*(angle->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 1))*(angle->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 2))*(angle->z);
    vprime3 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 0))*(angle->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 1))*(angle->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 2))*(angle->z);

    angle->x = vprime1;
    angle->y = vprime2;
    angle->z = vprime3;

    // position
    vprime1 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 0))*(position->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 1))*(position->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 2))*(position->z);
    vprime2 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 0))*(position->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 1))*(position->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 2))*(position->z);
    vprime3 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 0))*(position->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 1))*(position->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 2))*(position->z);

    position->x = vprime1;
    position->y = vprime2;
    position->z = vprime3;

    position->z = position->z + surface.height[surfaceIndex];
    position->x = position->x + surface.centerx[surfaceIndex];
    position->y = position->y + surface.centery[surfaceIndex];

}

// transform back to lab frame
void Image::transformInverse( Vector *position, Vector *angle, long surfaceIndex) {

    double vprime1, vprime2, vprime3;


    position->z = position->z - surface.height[surfaceIndex];
    position->x = position->x - surface.centerx[surfaceIndex];
    position->y = position->y - surface.centery[surfaceIndex];

    // angles
    vprime1 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 0))*(angle->x)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 1))*(angle->y)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 2))*(angle->z);
    vprime2 = (*(perturbation.inverserotationmatrix+9*surfaceIndex + 1*3 + 0))*(angle->x)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 1*3 + 1))*(angle->y)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 1*3 + 2))*(angle->z);
    vprime3 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 0))*(angle->x)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 1))*(angle->y)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 2))*(angle->z);

    angle->x = vprime1;
    angle->y = vprime2;
    angle->z = vprime3;

    // position
    vprime1 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 0))*(position->x) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 1))*(position->y) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 2))*(position->z);
    vprime2 = (*(perturbation.inverserotationmatrix+9*surfaceIndex + 1*3 + 0))*(position->x) +
        (*(perturbation.inverserotationmatrix+9*surfaceIndex + 1*3 + 1))*(position->y) +
        (*(perturbation.inverserotationmatrix+9*surfaceIndex + 1*3 + 2))*(position->z);
    vprime3 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 0))*(position->x) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 1))*(position->y) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 2))*(position->z);

    position->x = vprime1;
    position->y = vprime2;
    position->z = vprime3;


    // defocus
    position->z = position->z + perturbation.defocus[surfaceIndex] + surface.height[surfaceIndex];

    // decenter
    position->x = position->x + perturbation.decenterX[surfaceIndex] + surface.centerx[surfaceIndex];
    position->y = position->y + perturbation.decenterY[surfaceIndex] + surface.centery[surfaceIndex];

}

void Image::interceptDerivatives(Vector *normal, Vector position, long surfaceIndex)
{

    double rxd, udmin = 0.0, wdmin = 0.0;
    long umin = 0, rx, wmin = 0;
    double normal3, normal2;
    double phi, r;
    double dx, dy;
    double normalpr;

    dx = position.x - surface.centerx[surfaceIndex];
    dy = position.y - surface.centery[surfaceIndex];

    r = sqrt(dx*dx + dy*dy);
    rx = find_linear(&surface.radius[SURFACE_POINTS*surfaceIndex], SURFACE_POINTS, r, &rxd);
    if (perturbation.zernikeflag[surfaceIndex] == 1) {
        wmin=find_linear(perturbation.zernike_r_grid, PERTURBATION_POINTS, r/surface.rmax[surfaceIndex], &wdmin);
        phi=atan2(dy, dx);
        if (phi < 0) phi += 2*M_PI;
        umin = find_linear(perturbation.zernike_phi_grid, PERTURBATION_POINTS, phi, &udmin);
    }

    normal3 = interpolate_linear(&surface.normal[SURFACE_POINTS*surfaceIndex], rx, rxd);
    if (perturbation.zernikeflag[surfaceIndex]== 0) {
        normal->x = -normal3*dx/r;
        normal->y = -normal3*dy/r;
        normal->z = 1.0;
    } else {
        normal2 = 0;
        normalpr = normal3;
        normal3+=(interpolate_bilinear(perturbation.zernike_summed_nr_p+surfaceIndex*PERTURBATION_POINTS*PERTURBATION_POINTS, PERTURBATION_POINTS, umin, udmin, wmin, wdmin))/surface.rmax[surfaceIndex];
        normal2+=interpolate_bilinear(perturbation.zernike_summed_np_r+surfaceIndex*PERTURBATION_POINTS*PERTURBATION_POINTS, PERTURBATION_POINTS, umin, udmin, wmin, wdmin);
        normal->x = -normal3*dx/r + dy*normal2/(r*r);
        normal->y = -normal3*dy/r - dx*normal2/(r*r);
        normal->z = 1.0;
    }
    normalize(normal);

}

int Image::chooseSurface (long *newSurf, long *oldSurf) {

    if (photon.direction==1) {
    trynextsurface:;
        if ((*newSurf)<=(*oldSurf)+1) {
            (*newSurf)--;
            if ((*newSurf)==(*oldSurf)) (*newSurf)=(*oldSurf)-1;
            if ((*newSurf)<=-1) (*newSurf) = (*oldSurf)+2;
        } else {
            (*newSurf)++;
        }
        if ((*newSurf) > nsurf-1) {
            return(1);
        } else {
            if (ghost[(*newSurf)]==1) goto trynextsurface;
            return(0);
        }
    } else {
    trynextsurfaceback:;
        if ((*newSurf)>=(*oldSurf)-1) {
            (*newSurf)++;
            if ((*newSurf)==(*oldSurf)) (*newSurf) = (*oldSurf)+1;
            if ((*newSurf)>nsurf-1) (*newSurf)=(*oldSurf)-2;
        } else {
            (*newSurf)--;
        }
        if ((*newSurf) <=-1) {
            return(1);
        } else {
            if (ghost[(*newSurf)]==1) goto trynextsurfaceback;
            return(0);
        }
    }

}


void Image::atmosphereDiffraction (Vector *angle) {

    double radius, cr, cc, tf;
    fftw_plan pb;
    long i, j = 0, ix, jx;
    double dx, dy;
    double norm;
    double seeing;

    seeing = (totalseeing + 1e-6)*pow(1/cos(zenith), 0.6)*photon.wavelengthFactor;
    norm = sqrt(0.0229*pow(0.98*(1e-4*photon.wavelength)/(seeing*ARCSEC),-5./3.))*3.75e8;

    for (i = 0; i < SCREEN_SIZE; i++) {
        for (j = 0; j < SCREEN_SIZE; j++) {
            radius = sqrt((i - SCREEN_SIZE/2 + 0.5)*(i - SCREEN_SIZE/2 + 0.5)+
                          (j - SCREEN_SIZE/2 + 0.5)*(j - SCREEN_SIZE/2 + 0.5))*screen.fine_sizeperpixel;
            dx = (i - SCREEN_SIZE/2+0.5)*screen.fine_sizeperpixel*cos(photon.shiftedAngle) +
                (j - SCREEN_SIZE/2+0.5)*screen.fine_sizeperpixel*sin(photon.shiftedAngle);
            dy = -(i - SCREEN_SIZE/2 + 0.5)*screen.fine_sizeperpixel*sin(photon.shiftedAngle) +
                (j - SCREEN_SIZE/2 + 0.5)*screen.fine_sizeperpixel*cos(photon.shiftedAngle);
            if (radius > surface.innerRadius[0] && radius < surface.outerRadius[0]) {
                screen.inscreen[SCREEN_SIZE*i + j][0] = cos(screen.phasescreen[i*SCREEN_SIZE + j]*norm);
                screen.inscreen[SCREEN_SIZE*i + j][1] = sin(screen.phasescreen[i*SCREEN_SIZE + j]*norm);
                // for (k = 0; k<nspid; k++) {
                //     if (obstruction.spider_type[k] = =1) {
                //         if (fabs(dx-obstruction.spider_center[k]) < obstruction.spider_width[k]) {
                //             screen.inscreen[SCREEN_SIZE*i+j][0] = 0.0;
                //             screen.inscreen[SCREEN_SIZE*i+j][1] = 0.0;
                //         }
                //     }
                //     if (obstruction.spider_type[k]==2) {
                //         if (fabs(dy-obstruction.spider_center[k]) < obstruction.spider_width[k]) {
                //             screen.inscreen[SCREEN_SIZE*i+j][0] = 0.0;
                //             screen.inscreen[SCREEN_SIZE*i+j][1] = 0.0;
                //         }
                //     }
                // }
            } else {
                screen.inscreen[SCREEN_SIZE*i + j][0] = 0.0;
                screen.inscreen[SCREEN_SIZE*i + j][1] = 0.0;
            }
        }
    }


    pb = fftw_plan_dft_2d(SCREEN_SIZE,SCREEN_SIZE, screen.inscreen, screen.outscreen, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(pb);
    fftw_destroy_plan(pb);

    for (i = 0; i<SCREEN_SIZE; i++) {
        for (j = 0; j<SCREEN_SIZE; j++) {
            ix = i+SCREEN_SIZE/2;
            jx = j+SCREEN_SIZE/2;
            ix = ix % (SCREEN_SIZE);
            jx = jx % (SCREEN_SIZE);
            screen.focalscreen[ix*SCREEN_SIZE+jx] = pow(screen.outscreen[i*SCREEN_SIZE+j][0], 2)+
                pow(screen.outscreen[i*SCREEN_SIZE+j][1], 2);
        }
    }


    tf = 0.0;
    for (i = 0; i < SCREEN_SIZE; i++) for (j = 0; j < SCREEN_SIZE; j++) tf += screen.focalscreen[i*SCREEN_SIZE + j];
    cr = RngDouble()*tf;
    cc = 0.0;
    for (i = 0; i<SCREEN_SIZE; i++) {
        for (j = 0; j<SCREEN_SIZE; j++) {
            cc+=screen.focalscreen[i*SCREEN_SIZE+j];
            if (cr<cc) goto breakphoton;
        }
    }

 breakphoton:;

    angle->x = angle->x + (i-SCREEN_SIZE/2+RngDouble()-0.5)*photon.wavelength*1e-3/(SCREEN_SIZE*screen.fine_sizeperpixel);
    angle->y = angle->y + (j-SCREEN_SIZE/2+RngDouble()-0.5)*photon.wavelength*1e-3/(SCREEN_SIZE*screen.fine_sizeperpixel);
    angle->z = smallAnglePupilNormalize(angle->x, angle->y);

    photon.prtime = photon.time;

}


int Image::bloom (int saturatedFlag) {

    long newxpos, oldxpos, stepadd, iii;
    long startstep = 0;
    long location;

    location = chip.nampx*(photon.yPos - miny) + (photon.xPos - minx);

    if (photon.xPos >= chip.midpoint) stepadd = 1; else stepadd = 0;

    if (*(satupmap + location) >= 0 || *(satdownmap + location) >= 0) {

        oldxpos = photon.xPos;
        if (*(satupmap + location) < 0) startstep = photon.xPos - *(satdownmap + location);
        if (*(satdownmap + location) < 0) startstep = *(satupmap + location) - photon.xPos;
        if (*(satupmap + location) >= 0 && *(satdownmap + location) >= 0) {
            if (photon.xPos - *(satdownmap + location) < *(satupmap + location) - photon.xPos)
                startstep = photon.xPos - *(satdownmap + location); else
                startstep = *(satupmap + location) - photon.xPos;
        }

        for (iii = startstep; iii < chip.midpoint; iii++) {

            if (*(satupmap + location) >= 0) {
                newxpos = oldxpos + iii;
                if (newxpos < (chip.midpoint + stepadd*chip.midpoint)) {
                    if (*(chip.focal_plane + chip.nampx*(photon.yPos - miny) + (newxpos - minx)) < well_depth) {
                        photon.xPos = newxpos;
                        if (saturatedFlag == 1) *(satupmap + location) = newxpos;
                        return(0);
                    } else {
                        if (newxpos > *(satupmap+location) && saturatedFlag == 1) *(satupmap + location) = newxpos;
                    }
                } else {
                    if (saturatedFlag == 1) *(satupmap+location) = -1;
                }
            }

            if (*(satdownmap+location) >= 0) {
                newxpos = oldxpos - iii;
                if (newxpos >= stepadd*chip.midpoint) {
                    if (*(chip.focal_plane + chip.nampx*(photon.yPos-miny) + (newxpos - minx)) < well_depth) {
                        photon.xPos = newxpos;
                        if (saturatedFlag == 1) *(satdownmap + location) = newxpos;
                        return(0);
                    } else {
                        if (newxpos < *(satdownmap + location) && saturatedFlag == 1) *(satdownmap + location) = newxpos;
                    }
                } else {
                    if (saturatedFlag == 1) *(satdownmap + location) = -1;
                }
            }

        }
    }


    return(1);

}



void Image::saturate (long source, Vector* largeAngle)

{
    long location, origlocation;
    long leftover;
    long minrad;

    leftover = photon.sourceOver_m;
    origlocation = chip.nampx*(photon.yPos - miny)+(photon.xPos - minx);
    location = origlocation;

    if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny && photon.yPos <= maxy) {

    rebloom:;
        *(chip.focal_plane + location) += leftover;
        if (*(chip.focal_plane + location) > well_depth) {
            leftover = static_cast<long> (*(chip.focal_plane + location) - well_depth);
            *(chip.focal_plane + location) = (float)well_depth;
            if (blooming == 1) {
                if (bloom(1)) goto fullysat;
                location = chip.nampx*(photon.yPos - miny) + (photon.xPos - minx);
                goto rebloom;
            }
        }

    fullysat:;

        if (*(chip.focal_plane + origlocation) >= well_depth) {

            if (photon.ghostFlag == 0 && sources.spatialtype[source] != 4 && sources.spatialtype[source] != 1) {
                minrad = (long)fabs((largeAngle->y)/DEGREE*platescale/pixsize) - 4;
                if (minrad == (long)(fabs(photon.sourceSaturationRadius) + 1)) {
                    photon.sourceSaturationRadius = (double)minrad;
                }
            }
        }


    } else {

        if (photon.ghostFlag == 0 && sources.spatialtype[source] != 4 && sources.spatialtype[source] != 1 && sources.type[source] >= 4) {
            minrad = 0;
            double deltaX, deltaY;
            deltaX = (largeAngle->x)/DEGREE*platescale/pixsize;
            deltaY = (largeAngle->y)/DEGREE*platescale/pixsize;
            if (photon.yPos - deltaY > maxy && photon.yPos > maxy) if (-deltaY - (photon.yPos - maxy) > minrad) minrad = (long)(-deltaY - (photon.yPos - maxy));
            if (photon.yPos - deltaY < miny && photon.yPos < miny) if (deltaY - (photon.yPos - miny) > minrad) minrad = (long)(deltaY - (photon.yPos - miny));
            if (photon.xPos - deltaX > maxx && photon.xPos > maxx) if (-deltaX - (photon.xPos - maxx) > minrad) minrad = (long)(-deltaX - (photon.xPos - maxx));
            if (photon.xPos - deltaX < minx && photon.xPos < minx) if (deltaX - (photon.xPos - minx) > minrad) minrad = (long)(deltaX - (photon.xPos - minx));
            if (minrad > photon.sourceSaturationRadius) {
                photon.sourceSaturationRadius = (double)minrad;
            }
        }


    }

}

int Image::findSurface (Vector angle, Vector position, double *distance, long surfaceIndex) {


    double zmin;
    long miss;

    double a, b, c, l;
    double disc, l1, l2;
    double tol;
    double discrepancy;
    double zv;
    double radcurv;
    int iter = 0;

    if (surface.radiusCurvature[surfaceIndex] == 0.0) {

        l = (surface.height[surfaceIndex] - position.z)/angle.z;
        miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex);
        discrepancy = fabs(position.z + angle.z*l - zv);
        while (discrepancy > RAYTRACE_TOLERANCE) {
            l = l + (zv - (position.z + angle.z*l))/angle.z;
            miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex);
            discrepancy = fabs(position.z + angle.z*l - zv);
            iter++;
            if (iter > 10) goto iterjump;
        }
        zmin = discrepancy;
        *distance = l;
    } else {
        radcurv = -surface.radiusCurvature[surfaceIndex];
        a = (angle.x*angle.x + angle.y*angle.y + (1 + surface.conic[surfaceIndex])*angle.z*angle.z)*(-1.0)/radcurv;
        b = (-2.0)/radcurv*((1 + surface.conic[surfaceIndex])*angle.z*(position.z - surface.height[surfaceIndex])+
                            angle.x*position.x+angle.y*position.y) - 2*angle.z;
        c = (-1.0)/radcurv*((1 + surface.conic[surfaceIndex])*(position.z*position.z - 2*position.z*surface.height[surfaceIndex] +
                                                             surface.height[surfaceIndex]*surface.height[surfaceIndex]) +
                            position.x*position.x + position.y*position.y) - 2*(position.z - surface.height[surfaceIndex]);
        disc = b*b - 4*a*c;
        //printf("%g %i\n",disc,surfaceIndex);
        if (disc > 0) {
            l1 = ((-b + sqrt(disc))/2/a);
            l2 = ((-b - sqrt(disc))/2/a);
            l = 0.0;
            if ((l2 > 0) && (l1 > 0)) {
                if (l2 < l1) l = l2; else l = l1;
            }
            if ((l2 > 0) && (l1 < 0)) l = l2;
            if ((l2 < 0) && (l1 > 0)) l = l1;

            if (l == 0.0) {
                tol = 40000.0;
                miss = goldenBisectSurface(l - tol, l, l + tol, &zmin, angle, position, distance, surfaceIndex);
            } else {
                miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex);
                discrepancy = fabs(position.z + angle.z*l - zv);
                while (discrepancy > RAYTRACE_TOLERANCE) {
                    l = l + (zv - (position.z+angle.z*l))/angle.z;
                    miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex);
                    discrepancy = fabs(position.z + angle.z*l - zv);
                    iter++;
                    if (iter > 10) goto iterjump;
                }
                zmin = discrepancy;
                *distance = l;
            }
        } else {
        iterjump:;
            l = 0.0;
            tol = 40000.0;
            miss = goldenBisectSurface(l - tol, l, l + tol, &zmin, angle, position, distance, surfaceIndex);
        }
    }

    if (zmin > 10.0 || miss == 1) {
        return(1);
    } else {
        return(0);
    }

}


int Image::goldenBisectSurface(double a, double b, double c, double *z, Vector position, Vector angle, double *distance, long surfaceIndex) {

    double f0, f1, f2, f3;
    double x0, x1, x2, x3;
    double zv;
    long miss;

    x0 = a;
    x3 = c;
    if (fabs(c - b) > fabs(b - a)) {
        x1 = b;
        x2 = (b + HALFSQ5M1*(c - b));
    } else {
        x2 = b;
        x1 = (b - HALFSQ5M1*(b - a));
    }

    miss = getIntercept(position.x + angle.x*x1, position.y + angle.y*x1, &zv, surfaceIndex);
    f1 = fabs(position.z + angle.z*x1 - zv);
    miss = getIntercept(position.x + angle.x*x2, position.y + angle.y*x2, &zv, surfaceIndex);
    f2 = fabs(position.z + angle.z*x2 - zv);

    while (fabs(x2 - x1) > RAYTRACE_TOLERANCE) {

        if (f2 < f1) {
            x0 = x1;
            x1 = x2;
            x2 = HALF3MSQ5*x1 + HALFSQ5M1*x3;
            miss = getIntercept(position.x + angle.x*x2, position.y + angle.y*x2, &zv, surfaceIndex);
            f0 = f1;
            f1 = f2;
            f2 = fabs(position.z + angle.z*x2 - zv);
        } else {
            x3 = x2;
            x2 = x1;
            x1 = HALF3MSQ5*x2 + HALFSQ5M1*x0;
            miss = getIntercept(position.x + angle.x*x1, position.y + angle.y*x1, &zv, surfaceIndex);
            f3 = f2;
            f2 = f1;
            f1 = fabs(position.z + angle.z*x1 - zv);
        }
    }

    if (f1 < f2) {
        *z = f1;
        *distance = x1;
    } else {
        *z = f2;
        *distance = x2;
    }

    return(miss);

}

int Image::getIntercept(double x, double y, double *z, long surfaceIndex) {

    double r, phi;
    double uu, ww, tt;
    double dx, dy;
    long ttint;
    long lSurfaceIndex;

    dx = x - surface.centerx[surfaceIndex];
    dy = y - surface.centery[surfaceIndex];
    r = sqrt(dx*dx + dy*dy);
    lSurfaceIndex=SURFACE_POINTS*surfaceIndex;
    if ((r > surface.radius[lSurfaceIndex + SURFACE_POINTS - 1]) ||
        (r < surface.radius[lSurfaceIndex + 0])) return(1);

    ttint = find_linear(&surface.radius[lSurfaceIndex], SURFACE_POINTS, r, &tt);
    *z = interpolate_linear(&surface.profile[lSurfaceIndex], ttint, tt);

    find(&surface.radiusArea[lSurfaceIndex], SURFACE_POINTS, r, &(photon.vvint));

    photon.wwint = find_linear(perturbation.zernike_r_grid, PERTURBATION_POINTS, r/surface.rmax[surfaceIndex], &ww);
    phi = atan2(dy, dx);
    if (phi < 0) phi += 2*M_PI;
    photon.uuint = find_linear(perturbation.zernike_phi_grid, PERTURBATION_POINTS, phi, &uu);
    if (perturbation.zernikeflag[surfaceIndex] == 1) {
        *z += interpolate_bilinear(perturbation.zernike_summed + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex,
                                   PERTURBATION_POINTS, photon.uuint, uu, photon.wwint, ww);
    }

    return(0);

}

int Image::getDeltaIntercept(double x, double y, double *z, long surfaceIndex) {

    double r, phi;
    double uu, ww;
    long uuint, wwint;
    double dx, dy;
    int miss;

    miss = 0;

    dx = x - surface.centerx[surfaceIndex];
    dy = y - surface.centery[surfaceIndex];
    r = sqrt(dx*dx + dy*dy);


    find(&surface.radiusArea[SURFACE_POINTS*surfaceIndex], SURFACE_POINTS, r, &(photon.vvint));
    if ((r > surface.radius[SURFACE_POINTS*surfaceIndex + SURFACE_POINTS - 1]) ||
        (r < surface.radius[SURFACE_POINTS*surfaceIndex + 0])) miss=1;
    *z = 0.0;

    if (perturbation.zernikeflag[surfaceIndex] == 1) {
        wwint = find_linear(perturbation.zernike_r_grid, PERTURBATION_POINTS, r/surface.rmax[surfaceIndex], &ww);
        phi = atan2(dy, dx);
        if (phi < 0) phi += 2*M_PI;
        uuint = find_linear(perturbation.zernike_phi_grid, PERTURBATION_POINTS, phi, &uu);
        *z += interpolate_bilinear(perturbation.zernike_summed + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex,
                                   PERTURBATION_POINTS, uuint, uu, wwint, ww);
    }

    return(miss);
}


int Image::getWavelengthTime (double *wavelength, double *time, long source) {

    double tempf1;
    long index;
    double dustvalue = 1.0;

    // select wavelength
    tempf1 = RngDouble();
    find(sed_c + sed_ptr[sources.sedptr[source]], sed_n[sources.sedptr[source]], tempf1, &index);
    if (sources.type[source] == 1 || (sources.type[source] == 0 && domewave != 0.0)) {
        *wavelength = *(sed_w + sed_ptr[sources.sedptr[source]] + index);
    } else {
        *wavelength = interpolate(sed_w + sed_ptr[sources.sedptr[source]],sed_c + sed_ptr[sources.sedptr[source]], tempf1, index);
    }
    *wavelength = *wavelength/1000.;

    // dust at source
    if (sources.dusttypez[source] != 0) {
        if (sources.dusttypez[source] == 1) dustvalue = dust.ccm(*wavelength, sources.dustparz[source][0], sources.dustparz[source][1]);
        if (sources.dusttypez[source] == 2) dustvalue = dust.calzetti(*wavelength, sources.dustparz[source][0], sources.dustparz[source][1]);
        if (RngDouble() > dustvalue) return(1);
    }

    // redshift photon
    *wavelength = (*wavelength)*(1 + sources.redshift[source]);
    if (*wavelength < 0.3 || *wavelength >= 1.2) return(1);

    // dust at z=0
    if (sources.dusttype[source] != 0) {
        if (sources.dusttype[source] == 1) dustvalue = dust.ccm(*wavelength, sources.dustpar[source][0], sources.dustpar[source][1]);
        if (sources.dusttype[source] == 2) dustvalue = dust.calzetti(*wavelength, sources.dustpar[source][0], sources.dustpar[source][1]);
        if (RngDouble() > dustvalue) return(1);
    }

    // choose time
    *time = (double)((RngDouble())*exptime);

    return(0);

}

void Image::getAngle (Vector *angle, double time, long source) {

    if (sources.spatialtype[source] == MOVINGPOINT) {
        angle->x = sources.vx[source] + (sources.spatialpar[source][0]*cos(rotationangle)+
                                         sources.spatialpar[source][1]*(-sin(rotationangle)))*
            (time - exptime/2.0 + timeoffset)*ARCSEC;
        angle->y = sources.vy[source] + (sources.spatialpar[source][0]*sin(rotationangle)+
                                         sources.spatialpar[source][1]*cos(rotationangle))*
            (time - exptime/2.0 + timeoffset)*ARCSEC;
    } else {
        angle->x = sources.vx[source];
        angle->y = sources.vy[source];
    }
    angle->z = smallAnglePupilNormalize(angle->x, angle->y);

}

void Image::getDeltaAngle (Vector *angle, Vector *position, long source) {

    double dx=0., dy=0.;

    if (sources.spatialtype[source] != POINT && sources.spatialtype[source] != MOVINGPOINT) {
        if (sources.spatialtype[source] == SERSIC2D) {
            galaxy.sersic2d(sources.spatialpar[source][0], sources.spatialpar[source][1],
                            sources.spatialpar[source][2]*DEGREE, sources.spatialpar[source][3],
                            &dx, &dy);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        } else if (sources.spatialtype[source] == GAUSSIAN) {
            dx = random_gaussian()*ARCSEC*sources.spatialpar[source][0];
            dy = random_gaussian()*ARCSEC*sources.spatialpar[source][0];
        } else if (sources.spatialtype[source] == PINHOLE) {
            double r = sqrt(RngDouble()*sources.spatialpar[source][3]*sources.spatialpar[source][3]);
            double phi = RngDouble()*2*M_PI;
            double xt = position->x - r*cos(phi) - sources.spatialpar[source][0];
            double yt = position->y - r*sin(phi) - sources.spatialpar[source][1];
            double zt = position->z - surface.height[0] - sources.spatialpar[source][2];
            double finiteDistanceR = sqrt(xt*xt + yt*yt + zt*zt);
            dx = -xt/finiteDistanceR;
            dy = yt/finiteDistanceR;
        } else if (sources.spatialtype[source] == SERSIC) {
            galaxy.sersic(sources.spatialpar[source][0], sources.spatialpar[source][1],
                          sources.spatialpar[source][2], sources.spatialpar[source][3]*DEGREE,
                          sources.spatialpar[source][4]*DEGREE, sources.spatialpar[source][5],
                          &dx,&dy);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        } else if (sources.spatialtype[source] == IMAGE) {
            float imagepixel, localcumulative;
            long ii, jj, selectedsky;
            selectedsky = (long)sources.spatialpar[source][2];
            imagepixel = (RngDouble()*cumulative[selectedsky]);
            localcumulative = 0;
            for (ii = 0; ii < naxesb[selectedsky][0]; ii++) {
                if (cumulativex[selectedsky][ii] > imagepixel) {ii--; goto jumpa;}
            }
            ii--;
        jumpa:;
            localcumulative = cumulativex[selectedsky][ii];
            for (jj = 0; jj < naxesb[selectedsky][1]; jj++) {
                localcumulative += (*(tempptr[selectedsky] + (ii)*naxesb[selectedsky][1] + jj));
                if (localcumulative > imagepixel) {jj--; goto jumpb;}
            }
            jj--;
        jumpb:;
            dx=((jj - (naxesb[selectedsky][1]/2.0))*cos(sources.spatialpar[source][1]*DEGREE)+
                (-(ii - (naxesb[selectedsky][0]/2.0)))*sin(sources.spatialpar[source][1]*DEGREE))*
                sources.spatialpar[source][0]*ARCSEC;
            dy=-((jj - (naxesb[selectedsky][1]/2.0))*(-sin(sources.spatialpar[source][1]*DEGREE))+
                (-(ii - (naxesb[selectedsky][0]/2.0)))*cos(sources.spatialpar[source][1]*DEGREE))*
                sources.spatialpar[source][0]*ARCSEC;
        }
        if (sources.gamma1[source] != 0.0 || sources.gamma2[source] != 0.0) {
            double dxp = dx*(1 + sources.gamma1[source]-sources.kappa[source]) - dy*(sources.gamma2[source]);
            double dyp = dy*(1 - sources.gamma1[source]-sources.kappa[source]) - dx*(sources.gamma2[source]);
            dx = dxp;
            dy = dyp;
        }
        angle->x = angle->x - dx*cos(rotationangle) - dy*(sin(rotationangle));
        angle->y = angle->y - dx*(sin(rotationangle)) + dy*cos(rotationangle);
        angle->z = smallAnglePupilNormalize(angle->x, angle->y);
    }


}



int Image::photonSiliconPropagate(Vector *angle, Vector *position, double lambda, Vector normal, double dh, long waveSurfaceIndex) {

    double travel,dtravel;
    long yindex;
    double ryindex;
    Vector origAngle;
    double dead;

    photon.z0 = position->z;
    vectorCopy(*angle, &origAngle);

    // silicon refraction
    yindex = find_linear(silicon.wavelengthGrid, silicon.numWavelength, photon.wavelength, &ryindex);
    double nSi = interpolate_linear(silicon.indexRefraction, yindex, ryindex);
    refract(angle, normal, 1, nSi);

    // photo-electron conversion
    photon.xindex = find_linear(silicon.temperatureGrid, silicon.numTemperature, ccdtemp, &(photon.rxindex));
    double mfp = interpolate_bilinear(silicon.meanFreePath, silicon.numWavelength, photon.xindex, photon.rxindex, yindex, ryindex);
    double randNum;
    double conversion;
    if (photoelectric == 1) {
        conversion = 1.0 - exp(-2*siliconthickness/1e3/fabs(angle->z)/mfp);
    } else {
        conversion = 1.0;
        mfp = 0.0;
    }
    (photon.counter)++;
    if (conversion > dynamicTransmission[2*natmospherefile + 2*nsurf + 1 + waveSurfaceIndex]) {
        dynamicTransmission[2*natmospherefile + 2*nsurf + 1 + waveSurfaceIndex] = conversion;
    }
    if (photon.counter <= photon.maxcounter) {
        randNum = saveRand[photon.counter];
    } else {
        randNum = RngDouble();
    }
    if (randNum > conversion) return(1);
    travel = mfp*(-log(1.0 - randNum));
    photon.location = chip.nampx*(photon.yPos - miny) + (photon.xPos - minx);
    if (photon.xPos < minx || photon.xPos > maxx || photon.yPos < miny || photon.yPos > maxy) {
        long xL, yL;
        xL = photon.xPos;
        yL = photon.yPos;
        if (photon.xPos < minx) xL = minx;
        if (photon.xPos > maxx) xL = maxx;
        if (photon.yPos < miny) yL = miny;
        if (photon.yPos > maxy) yL = maxy;
        photon.location = chip.nampx*(yL - miny) + (xL - minx);
    }
    if (deadlayer == 1) {
        dead = silicon.deadLayer[photon.location];
    } else {
        dead = 0.0;
    }
    if (fabs(travel*angle->z) < dead/1e6) return(1);
    photon.collect_z = photon.z0 + angle->z/fabs(angle->z)*siliconthickness/1e3;
    if (fabs(travel*angle->z) >= siliconthickness/1e3) {
        double reflectivity = fringing(origAngle, normal, photon.wavelength, nSi, (double)siliconthickness + dh*1000.0);
        if (RngDouble() > reflectivity || fringeflag == 0) return(1);
        photon.ghostFlag = 1;
        dtravel = fabs(siliconthickness/1e3/angle->z);
        propagate(position, *angle, dtravel);
        reflect(angle, normal);
        propagate(position, *angle, travel - dtravel);
    } else {
        propagate(position, *angle, travel);
    }
    return(0);
}

int Image::electronSiliconPropagate(Vector *angle, Vector *position) {

    long windex, zindex, uindex;
    double rwindex, rzindex, ruindex;
    double dopant;

    // charge diffusion
    zindex = find_linear(silicon.thicknessGrid, silicon.numThickness, siliconthickness/1e4 - fabs((position->z - photon.z0)/10.0), &rzindex);
    if (impurityvariation == 1) {
        dopant = silicon.nbulkmap[photon.location]*nbulk;
    } else {
        dopant = nbulk;
    }
    windex = find_linear(silicon.dopantGrid, silicon.numDopant, dopant, &rwindex);
    double sg = interpolate_trilinear(silicon.sigma, silicon.numTemperature, silicon.numThickness,
                                      windex, rwindex, photon.xindex, photon.rxindex, zindex, rzindex);
    if (chargediffusion == 0) sg = 0.0;

    // complications to charge diffusion (effect of lateral fields)
    double fsg = interpolate_bilinear_float(silicon.fsigma, silicon.numThickness, windex, rwindex, zindex, rzindex);
    double gsg = interpolate_bilinear_float(silicon.gsigma, silicon.numThickness, windex, rwindex, zindex, rzindex);
    double sa, sb, ga, gb, da, db;
    double chsx, chsy, chv, chp;
    sa = silicon.sigmaX[photon.location];
    sb = silicon.sigmaY[photon.location];
    ga = silicon.gammaX[photon.location];
    gb = silicon.gammaY[photon.location];
    da = silicon.deltaX[photon.location];
    db = silicon.deltaY[photon.location];
    if (chargesharing == 1) {
        chsx = 0.0;
        chsy = 0.0;
        double rho;
        if (photon.xPos <= maxx && photon.xPos >= minx && photon.yPos >= miny && photon.yPos <= maxy) {
            rho = sqrt((photon.xPosR)*(photon.xPosR) + (photon.yPosR)*(photon.yPosR));
            uindex = find_linear(silicon.rho, silicon.numTemperature, rho*pixsize*1e-4, &ruindex);
            chp = *(chip.focal_plane + chip.nampx*(photon.yPos - miny) + (photon.xPos - minx));
            chv = interpolate_trilinear(silicon.hsigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsx += chp*chv*(photon.xPosR)/rho;
            chsy += chp*chv*(photon.yPosR)/rho;
        }
        if (photon.xPos + 1 <= maxx && photon.xPos + 1 >= minx && photon.yPos >= miny && photon.yPos <= maxy) {
            rho = sqrt((photon.xPosR - 1)*(photon.xPosR - 1) + (photon.yPosR)*(photon.yPosR));
            uindex = find_linear(silicon.rho, silicon.numTemperature, rho*pixsize*1e-4, &ruindex);
            chp = *(chip.focal_plane + chip.nampx*(photon.yPos - miny) + (photon.xPos + 1 - minx));
            chv = interpolate_trilinear(silicon.hsigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsx += chp*chv*(photon.xPosR - 1)/rho;
            chsy += chp*chv*(photon.yPosR)/rho;
        }
        if (photon.xPos - 1 <= maxx && photon.xPos - 1 >= minx && photon.yPos >= miny && photon.yPos <= maxy) {
            rho = sqrt((photon.xPosR + 1)*(photon.xPosR + 1)+(photon.yPosR)*(photon.yPosR));
            uindex = find_linear(silicon.rho, silicon.numTemperature, rho*pixsize*1e-4, &ruindex);
            chp = *(chip.focal_plane + chip.nampx*(photon.yPos - miny) + (photon.xPos - 1 - minx));
            chv = interpolate_trilinear(silicon.hsigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsx += chp*chv*(photon.xPosR + 1)/rho;
            chsy += chp*chv*(photon.yPosR)/rho;
        }
        if (photon.xPos <= maxx && photon.xPos >= minx && photon.yPos + 1 >= miny && photon.yPos + 1 <= maxy) {
            rho = sqrt((photon.xPosR)*(photon.xPosR)+(photon.yPosR - 1)*(photon.yPosR - 1));
            uindex = find_linear(silicon.rho, silicon.numTemperature, rho*pixsize*1e-4, &ruindex);
            chp = *(chip.focal_plane + chip.nampx*(photon.yPos + 1 - miny) + (photon.xPos - minx));
            chv = interpolate_trilinear(silicon.hsigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsx += chp*chv*(photon.xPosR)/rho;
            chsy += chp*chv*(photon.yPosR -1)/rho;
        }
        if (photon.xPos <= maxx && photon.xPos >= minx && photon.yPos - 1 >= miny && photon.yPos - 1 <= maxy) {
            rho = sqrt((photon.xPosR)*(photon.xPosR)+(photon.yPosR + 1)*(photon.yPosR + 1));
            uindex = find_linear(silicon.rho, silicon.numTemperature, rho*pixsize*1e-4, &ruindex);
            chp = *(chip.focal_plane + chip.nampx*(photon.yPos - 1 - miny) + (photon.xPos - minx));
            chv = interpolate_trilinear(silicon.hsigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsx += chp*chv*(photon.xPosR)/rho;
            chsy += chp*chv*(photon.yPosR + 1)/rho;
        }
        if (photon.xPos <= maxx && photon.xPos >= minx && photon.yPos >= miny && photon.yPos - 1 <= maxy) {
            rho = fabs(photon.yPosR + 0.5);
            uindex = find_linear(silicon.rho, silicon.numTemperature, rho*pixsize*1e-4, &ruindex);
            chp = silicon.chargeStopCharge;
            chv = interpolate_trilinear(silicon.isigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsy += chp*chv*(photon.yPosR + 0.5)/rho;
        }
        if (photon.xPos <= maxx && photon.xPos >= minx && photon.yPos + 1 >= miny && photon.yPos <= maxy) {
            rho = fabs(photon.yPosR - 0.5);
            uindex = find_linear(silicon.rho, silicon.numTemperature, rho*pixsize*1e-4, &ruindex);
            chp = silicon.chargeStopCharge;
            chv = interpolate_trilinear(silicon.isigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsy += chp*chv*(photon.yPosR -0.5)/rho;
        }
    } else {
        chsx = 0.0;
        chsy = 0.0;
    }
    if (fieldanisotropy == 0) {
        sa = 0.0;
        sb = 0.0;
        da = 0.0;
        db = 0.0;
    }
    if (pixelerror == 0) {
        ga = 0.0;
        gb = 0.0;
    }
    position->x += (sg*random_gaussian() + fsg*sa + gsg*da + ga + chsx)/1e3;
    position->y += (sg*random_gaussian() + fsg*sb + gsg*db + gb + chsy)/1e3;
    position->z = photon.collect_z;
    return(0);

}

double Image::fringing (Vector angle, Vector normal, double wavelength, double nSi, double thickness) {

    double arg = fabs(normal.x*angle.x + normal.y*angle.y + normal.z*angle.z);
    if (arg > 1) arg = 1.0;
    double airAngle = acos(arg);
    int polarization = -1;
    if (RngDouble() < 0.5) polarization = 1;
    double siliconAngle = asin(sin(airAngle)/nSi);
    double airAngleOut = -airAngle;
    double delta = 2*M_PI/wavelength*thickness*nSi*sqrt(1.0 - sin(airAngle)*sin(airAngle)/(nSi*nSi));
    std::complex<double> rhot1 ((pow(cos(airAngle), polarization) - nSi*pow(cos(siliconAngle), polarization))/
                                (pow(cos(airAngle), polarization) + nSi*pow(cos(siliconAngle), polarization)), 0.0);
    std::complex<double> rhot2 ((nSi*pow(cos(siliconAngle), polarization) - pow(cos(airAngleOut), polarization))/
                                (pow(cos(airAngleOut), polarization) + nSi*pow(cos(siliconAngle), polarization)), 0.0);
    std::complex<double> phase (cos(-2.0*delta), sin(-2.0*delta));
    std::complex<double> gamma;
    std::complex<double> one (1, 0);
    gamma = (rhot1 + rhot2*phase)/(one + rhot1*rhot2*phase);
    double reflection = real(gamma*conj(gamma));
    return(reflection);


}


int Image::contaminationSurfaceCheck(Vector position, Vector *angle, long surfaceIndex) {

    if (RngDouble() > (double)(*(contamination.transmission +
                               PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex +
                               photon.vvint*PERTURBATION_POINTS + photon.uuint))) {
        return(1);
    }
    if (*(contamination.surfacelistmap + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex +
          photon.vvint*PERTURBATION_POINTS + photon.uuint) != -1) {
        long cc = *(contamination.surfacelistmap + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex +
                    photon.vvint*PERTURBATION_POINTS + photon.uuint);

        if (sqrt(pow(position.x - contamination.surfacelistx[surfaceIndex][cc], 2.0)+
                 pow(position.y - contamination.surfacelisty[surfaceIndex][cc], 2.0)) <
            contamination.surfacelists[surfaceIndex][cc]) {

            if (RngDouble() > exp(-contamination.absorptionLength*contamination.surfacelists[surfaceIndex][cc])) return(1);
            long index;
            double mu;
            if (contamination.surfacelistt[surfaceIndex][cc] == 0) {
                find(contamination.henyey_greenstein, contamination.elements, RngDouble(), &index);
                mu = contamination.henyey_greenstein_mu[index];
            } else {
                find(contamination.henyey_greenstein_w, contamination.elements, RngDouble(), &index);
                mu = contamination.henyey_greenstein_mu_w[index];
            }
            double phi = RngDouble()*2.0*M_PI;
            shift_mu(angle, mu, phi);
            if (mu < 0) {
                return(1);
            }
        }
    }
    return(0);


}
