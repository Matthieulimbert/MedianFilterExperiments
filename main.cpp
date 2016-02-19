//
//  main.cpp
//  MedianFilterExperiments
//
//  Created by Matthieu Limbert on 27/12/15.
//  Copyright © 2015 Matthieu Limbert. All rights reserved.
//

#include <iostream>
#include <opencv2/opencv.hpp>
#include "functions.hpp"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

int main(int argc, const char * argv[]) {

    
    // Definition of some parameters
    int n = 200;
    double p = 0.2;
    int rMax = 12;
    int numberOfIterations = 1;
    int numberOfOrientations = 8;
    
    
    // we create the noiseless image
    cv::Mat noiselessImage(n, n, cv::DataType<bool>::type);
    createNoiselessImage(noiselessImage);
    
    // we show the noiseless image
    showImage("noiseless image", noiselessImage);
    
    
    for (int iterationIndex = 0; iterationIndex < numberOfIterations; iterationIndex++) {
    
        // we create the noisy image
        cv::Mat noisyImage(n, n, cv::DataType<bool>::type);
        createNoisyImage(noiselessImage, noisyImage, p);
        
        
        // we show the noisy image
        showImage("noisy image", noisyImage);
        
        
        // we estimate the Bernoulli noise
        float pEstimation = noiseEstimation(noisyImage, rMax);
        
        
        // we estimate the contamination rate and the distance to the discontinuity
        cv::Mat contaminationMap(n, n, cv::DataType<float>::type);
        cv::Mat distanceMap(n, n, cv::DataType<int>::type);
        computeContaminationAndDistanceMaps(noisyImage, rMax, pEstimation, contaminationMap, distanceMap);
        
        
        // we show the contamination map
        showContaminationMap("contamination map", contaminationMap);
        
        
        // we denoise the image with a fixed circle window size median filter
        cv::Mat denoisedImage1(n, n, cv::DataType<bool>::type);
        denoiseFixedCircleMF(noisyImage, rMax, denoisedImage1);
        //showImage("denoised image with a fixed circle window median filter", denoisedImage1);
        showImage("denoisedImage1", denoisedImage1);
        
        
        // we denoise the image with an adaptive circle window size median filter
        cv::Mat denoisedImage2(n, n, cv::DataType<bool>::type);
        denoiseAdaptiveCircleMF(noisyImage, distanceMap, rMax, denoisedImage2);
        //showImage("denoised image with an adaptive circle window median filter", denoisedImage2);
        showImage("denoisedImage2", denoisedImage2);
        
        
        // we denoise the image with an adaptive elliptic window size median filter
        cv::Mat denoisedImage3(n, n, cv::DataType<bool>::type);
        denoiseAdaptiveEllipticMF(noisyImage, distanceMap, pEstimation, rMax, numberOfOrientations, denoisedImage3);
        //showImage("denoised image with an adaptive ellipitic window median filter", denoisedImage3);
        showImage("denoisedImage3", denoisedImage3);
        
        
        // we compute the error rates
        float errorRate1 = computeErrorRate(noiselessImage, denoisedImage1);
        float errorRate2 = computeErrorRate(noiselessImage, denoisedImage2);
        float errorRate3 = computeErrorRate(noiselessImage, denoisedImage3);
        
        // we update the means and the medians of the error rates
    
    }
    
    
    // we write the means and the medians of the error rate for the current Bernoulli noise
    
    
    

    return 0;
}
