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
    float pMin = 0.01f;
    float pMax = 0.49f;
    float pStep = 0.01f;
    int rMax = 12;
    int numberOfIterations = 1000;
    int numberOfOrientations = 8;
    
    
    // we create the noiseless image
    cv::Mat noiselessImage(n, n, cv::DataType<bool>::type);
    createNoiselessImage(noiselessImage);
    
    // we show the noiseless image
    //showImage("noiseless image", noiselessImage);
    
    // initialization of the vectors containing the results
    std::vector<float> means1, means2, means3;
    std::vector<float> medians1, medians2, medians3;
    std::vector<float> noiseLevels;
    
    std::cout << "Start median filters experiments... " << std::endl;
    
    // loop on Bernoulli noise level
    for (float p = pMin; p <= pMax; p += pStep ) {
        
        // display some information
        std::cout << std::endl << "noise level : " << p << std::endl;
        
        // update vector containing the noise levels
        noiseLevels.push_back(p);
        
        // Loop on the iterations
        std::vector<float> currentNoiseErrorRateVector1, currentNoiseErrorRateVector2, currentNoiseErrorRateVector3;
        for (int iterationIndex = 0; iterationIndex < numberOfIterations; iterationIndex++) {
            
            // display some information
            std::cout << "iteration : " << iterationIndex << std::endl;
            
            // we create the noisy image
            cv::Mat noisyImage(n, n, cv::DataType<bool>::type);
            createNoisyImage(noiselessImage, noisyImage, p);
            
            
            // we show the noisy image
            //showImage("noisy image", noisyImage);
            
            
            // we estimate the Bernoulli noise
            float pEstimation = noiseEstimation(noisyImage, rMax);
            
            
            // we estimate the contamination rate and the distance to the discontinuity
            cv::Mat contaminationMap(n, n, cv::DataType<float>::type);
            cv::Mat distanceMap(n, n, cv::DataType<int>::type);
            computeContaminationAndDistanceMaps(noisyImage, rMax, pEstimation, contaminationMap, distanceMap);
            
            
            // we show the contamination map
            //showContaminationMap("contamination map", contaminationMap);
            
            
            // we denoise the image with a fixed circle window size median filter
            cv::Mat denoisedImage1(n, n, cv::DataType<bool>::type);
            denoiseFixedCircleMF(noisyImage, rMax, denoisedImage1);
            //showImage("denoisedImage1", denoisedImage1);
            
            
            // we denoise the image with an adaptive circle window size median filter
            cv::Mat denoisedImage2(n, n, cv::DataType<bool>::type);
            denoiseAdaptiveCircleMF(noisyImage, distanceMap, rMax, denoisedImage2);
            //showImage("denoisedImage2", denoisedImage2);
            
            
            // we denoise the image with an adaptive elliptic window size median filter
            cv::Mat denoisedImage3(n, n, cv::DataType<bool>::type);
            denoiseAdaptiveEllipticMF(noisyImage, distanceMap, pEstimation, rMax, numberOfOrientations, denoisedImage3);
            //showImage("denoisedImage3", denoisedImage3);
            
            
            // we compute the error rates
            float errorRate1 = computeErrorRate(noiselessImage, denoisedImage1);
            float errorRate2 = computeErrorRate(noiselessImage, denoisedImage2);
            float errorRate3 = computeErrorRate(noiselessImage, denoisedImage3);
            currentNoiseErrorRateVector1.push_back(errorRate1);
            currentNoiseErrorRateVector2.push_back(errorRate2);
            currentNoiseErrorRateVector3.push_back(errorRate3);
            
            
        }
        
        // we write the means of the error rate for the current Bernoulli noise
        float currentNoiseMean1 = computeMean(currentNoiseErrorRateVector1);
        float currentNoiseMean2 = computeMean(currentNoiseErrorRateVector2);
        float currentNoiseMean3 = computeMean(currentNoiseErrorRateVector3);
        means1.push_back(currentNoiseMean1);
        means2.push_back(currentNoiseMean2);
        means3.push_back(currentNoiseMean3);
        
        // we write the medians of the error rate for the current Bernoulli noise
        float currentNoiseMedian1 = computeMedian(currentNoiseErrorRateVector1);
        float currentNoiseMedian2 = computeMedian(currentNoiseErrorRateVector2);
        float currentNoiseMedian3 = computeMedian(currentNoiseErrorRateVector3);
        medians1.push_back(currentNoiseMedian1);
        medians2.push_back(currentNoiseMedian2);
        medians3.push_back(currentNoiseMedian3);
        
        
    }
    
    // store the results in an excel sheet
    // TO DO
    std::string filename = "/Users/mlimbert/Documents/MedianFilterExperiments/ResultsNoiseLevels.xml";
    cv::FileStorage fs(filename, cv::FileStorage::WRITE);
    //fs.open(filename, cv::FileStorage::WRITE);
    fs << "noiseLevels" << noiseLevels;
    fs.release();
    
    filename = "/Users/mlimbert/Documents/MedianFilterExperiments/ResultsMeans1.xml";
    fs.open(filename, cv::FileStorage::WRITE);
    fs << "means1" << means1;
    fs.release();
    
    filename = "/Users/mlimbert/Documents/MedianFilterExperiments/ResultsMeans2.xml";
    fs.open(filename, cv::FileStorage::WRITE);
    fs << "means2" << means2;
    fs.release();
    
    filename = "/Users/mlimbert/Documents/MedianFilterExperiments/ResultsMeans3.xml";
    fs.open(filename, cv::FileStorage::WRITE);
    fs << "means3" << means3;
    fs.release();
    
    filename = "/Users/mlimbert/Documents/MedianFilterExperiments/ResultsMedians1.xml";
    fs.open(filename, cv::FileStorage::WRITE);
    fs << "medians1" << medians1;
    fs.release();
    
    filename = "/Users/mlimbert/Documents/MedianFilterExperiments/ResultsMedians2.xml";
    fs.open(filename, cv::FileStorage::WRITE);
    fs << "medians2" << medians2;
    fs.release();
    
    filename = "/Users/mlimbert/Documents/MedianFilterExperiments/ResultsMedians3.xml";
    fs.open(filename, cv::FileStorage::WRITE);
    fs << "medians3" << medians3;
    fs.release();
    

    
    
    
//        fs << "means1" << means1;
//        fs << "means2" << means2;
//        fs << "means3" << means3;
//        fs << "medians1" << medians1;
//        fs << "medians2" << medians2;
//        fs << "medians3" << medians3;

    
    
    return 0;
}
