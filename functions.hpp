//
//  functions.hpp
//  MedianFilterExperiments
//
//  Created by Matthieu Limbert on 29/12/15.
//  Copyright © 2015 Matthieu Limbert. All rights reserved.
//

#ifndef functions_hpp
#define functions_hpp
#include <opencv2/opencv.hpp>

#include <stdio.h>

// functions to create the noiseless image and the noisy images
void createNoiselessImage(cv::Mat &imageInOut);
void createNoisyImage(cv::Mat &noiselessImageIn, cv::Mat &noisyImageOut, float noise);
bool insideEllipseTest(int i, int j, int sizeImageIn);

// function to display an image
void showImage(std::string windowName, cv::Mat binaryImageIn);
void showContaminationMap(std::string windowName, cv::Mat contaminationMapIn);

// function to estimate the noise
float noiseEstimation(cv::Mat &noisyImage, int rMax);

// function that estimates the
void computeContaminationAndDistanceMaps(cv::Mat &noisyImageIn, int rMaxIn, float pEstimation, cv::Mat &contaminationMapOut, cv::Mat &distanceMapOut);

// Denoising functions
void denoiseFixedCircleMF(cv::Mat &noisyImageIn, int rMaxIn, cv::Mat &denoisedImageOut);
void denoiseAdaptiveCircleMF(cv::Mat const &noisyImageIn, cv::Mat const &distanceMapIn, int rMaxIn, cv::Mat &denoisedImageOut);
void denoiseAdaptiveEllipticMF(cv::Mat const &noisyImageIn, cv::Mat const &distanceMapIn, float pEstimationIn, int rMaxIn, int numberOfOrientationsIn, cv::Mat &denoisedImageOut);

// computation of the error rate
float computeErrorRate(cv::Mat const &image1In, cv::Mat const &image2In);

#endif /* functions_hpp */
