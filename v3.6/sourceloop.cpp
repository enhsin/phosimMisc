///
/// @package phosim
/// @file sourceloop.cpp
/// @brief source loop
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///
#include <pthread.h>

int Image::sourceLoop() {

    char tempstring[4096];
    long sourceCounter;
    long surfaceLimit;
    pthread_t *thread;
    thread_args *args;

    //    OPTIMIZATION INITIALIZATION
    thread = (pthread_t*)malloc(numthread*sizeof(pthread_t));
    args = (thread_args*)malloc(numthread*sizeof(thread_args));
    surfaceLimit = natmospherefile*2 + nsurf*2 + 2;
    state.dynamicTransmission = static_cast<double*>(malloc(surfaceLimit*(MAX_WAVELENGTH-MIN_WAVELENGTH+1)*sizeof(double)));
    for (int i = 0; i < MAX_WAVELENGTH-MIN_WAVELENGTH+1; i++) {
        for (int j = 0; j < surfaceLimit; j++) {
            state.dynamicTransmission[i*surfaceLimit + j] = -1.0;
        }
    }
    if (checkpointcount != 0) readCheckpoint(checkpointcount);

    //    LOGGING INITIALIZATION
    if (eventfile) {
        eventFitsFileName = std::string(this->outputfilename);
        eventFitsFileName = eventFitsFileName.replace(eventFitsFileName.find("_e_"), 3, "_r_") + ".fits";
        state.pEventLogging = new EventFile((int)(nphot*100), outputdir, eventFitsFileName);
    } else {
        state.pEventLogging = NULL;
    }
    pGrating = new Grating();
    if (throughputfile) initThroughput(&state.throughputLog, nsurf);
    counterInit(&state.counterLog);
    counterClear(&state.globalLog);
    state.cx = 0.0;
    state.cy = 0.0;
    state.cz = 0.0;
    state.r0 = 0.0;
    state.epR = 0.0;
    pthread_mutex_init(&lock.lock1, NULL);
    pthread_mutex_init(&lock.lock2, NULL);
    pthread_mutex_init(&lock.lock3, NULL);


    //    SOURCE TYPE LOOP
    for (int sourceType = 0; sourceType < 216; sourceType++) {

        sourceCounter = 0;

        if (static_cast<int>(sourceType/216.0*(checkpointtotal+1)) == checkpointcount) {

            //    SOURCE LOOP
            long subsource = 0;

            for (long source = 0; source < nsource; source++) {

                        if ((sourceType <= 214 && sources.type[source] == 0 && (source%215) == sourceType) ||
                            (sourceType == 215 && sources.type[source] > 0)) {

                            args[subsource].instance = this;
                            args[subsource].ssource = source;
                            pthread_create(&thread[subsource], NULL, &Image::threadFunction, &args[subsource]);
                            subsource++;
                            sourceCounter++;

                        }
                        if (subsource >= numthread || source == nsource - 1) {
                            for (long ss = 0; ss < subsource; ss++) {
                                pthread_join(thread[ss], NULL);
                            }
                            subsource = 0;
                        }

            }

        }
        if (sourceType >= 0  && sourceType < 215) sprintf(tempstring, "Dome Light         ");
        if (sourceType == 215) sprintf(tempstring, "Others ");
        if (sourceCounter > 0) counterCheck(&state.counterLog, sourceCounter, tempstring);
    }

    // COSMIC RAYS
    long long ray;
    long long detRay;
    if (checkpointcount == checkpointtotal) {
        if (backgroundMode > 0) {
            detRay = 0;
            ray = 0;
            cosmicRays(&detRay);
            if (detRay > 0) {
                for (long i = 0; i < detRay; i++) {
                    countGood(&state.counterLog, 1, &ray);
                }
                sprintf(tempstring, "Cosmic Rays        ");
                counterCheck(&state.counterLog, detRay, tempstring);
            }
        }
    }

    // OUTPUT DATA
    if (checkpointcount == checkpointtotal) writeImageFile();
    if (opdfile) writeOPD();
    if (checkpointcount !=  checkpointtotal) writeCheckpoint(checkpointcount);
    if (centroidfile) writeCentroidFile(outputdir, outputfilename, sourcePhoton, sourceXpos, sourceYpos, sources.id, nsource);
    if (throughputfile) writeThroughputFile(outputdir, outputfilename, &state.throughputLog, nsurf);
    if (eventfile) state.pEventLogging->eventFileClose();
    return(0);


}
