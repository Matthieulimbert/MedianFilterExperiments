//
//  functions.cpp
//  MedianFilterExperiments
//
//  Created by Matthieu Limbert on 29/12/15.
//  Copyright © 2015 Matthieu Limbert. All rights reserved.
//

#include "functions.hpp"


void createNoisyImage(cv::Mat &noiselessImageIn, cv::Mat &noisyImageOut, float noise)
{
    /*
     This function creates a noiseless white image containing a black ellipse on the foreground
     */
    
    assert(noiselessImageIn.rows == noiselessImageIn.cols); // noiseless image is a square
    assert(noisyImageOut.rows == noisyImageOut.cols);       // noisy image is a square
    assert(noiselessImageIn.rows == noisyImageOut.rows);    // noiseless and noisy images are same size
    int n = noiselessImageIn.rows;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
        {
            int a = rand();
            float b = (float)a / RAND_MAX;
            
            if (b >= noise)
            {
                noisyImageOut.at<bool>(i,j) = noiselessImageIn.at<bool>(i,j);
            }
            
            else if (b < noise && noiselessImageIn.at<bool>(i,j) == true)
            {
                noisyImageOut.at<bool>(i,j) = false;
            }
            
            else if (b < noise && noiselessImageIn.at<bool>(i,j) == false)
            {
                noisyImageOut.at<bool>(i,j) = true;
            }
            
        }
        
    }
}


void createNoiselessImage(cv::Mat &imageInOut)
{
    /*
    This function creates a noiseless white image containing a black ellipse on the foreground
    */
    
    assert(imageInOut.rows == imageInOut.cols);
    int n = imageInOut.rows;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
        {
            bool insideEllipse = insideEllipseTest(i, j, n);
            if (insideEllipse == true) {
                imageInOut.at<bool>(i,j) = true; // The pixel is inside the ellipse (black)
            }
            
            else {
                imageInOut.at<bool>(i,j) = false; // The pixel is outside the ellipse (white)
            }
            
        }
    }
}


bool insideEllipseTest(int i, int j, int sizeImageIn)
{
    // Initialisation of variables :
    float ellipse_1;
    float ellipse_2;
    float ellipse_3;
    float ellipse;
    bool testResult;
    
    // Calculation of whether the pixel is inside or outside of the ellipse :
    ellipse_1 = (float)sizeImageIn / 2 - j - 1;
    ellipse_2 = (float)sizeImageIn / 7;
    ellipse_1 = ellipse_1 / ellipse_2;
    ellipse_1 = ellipse_1 * ellipse_1;
    
    ellipse_2 = (float)sizeImageIn / 2 - i - 1;
    ellipse_3 = (float)sizeImageIn / 5;
    ellipse_2 = ellipse_2 / ellipse_3;
    ellipse_2 = ellipse_2 * ellipse_2;
    
    ellipse = ellipse_1 + ellipse_2;
    if (ellipse <= 1) {
        testResult = true; // The pixel is inside the ellipse
    }
    else {
        testResult = false; // The pixel is outside the ellipse
    }
    
    
    return testResult;
}


void showImage(std::string windowName, cv::Mat binaryImageIn)
{
    // we create the grayscale image to display it later
    cv::Mat grayscaleImage(binaryImageIn.size(), cv::DataType<uchar>::type);
    
    // we fill in the grayscale image
    for (int i = 0; i < binaryImageIn.rows; i++) {
        for (int j = 0; j < binaryImageIn.cols; j++)
        {
            if(binaryImageIn.at<bool>(i,j) == true)
            {
                grayscaleImage.at<uchar>(i,j) = 0; // The pixel is black
            }
            
            else
            {
                grayscaleImage.at<uchar>(i,j) = 255; // The pixel is white
            }
        }
    }
    
   // we display the grayscale image
    cv::namedWindow(windowName, cv::WINDOW_AUTOSIZE );
    imshow(windowName, grayscaleImage);
}


void showContaminationMap(std::string windowName, cv::Mat contaminationMapIn)
{
    // we create the grayscale image to display it later
    cv::Mat grayscaleContaminationMap(contaminationMapIn.size(), cv::DataType<uchar>::type);
    
    // we fill in the grayscale image
    for (int i = 0; i < contaminationMapIn.rows; i++) {
        for (int j = 0; j < contaminationMapIn.cols; j++)
        {
            // we make sure the contamination map has a the right structure
            assert(contaminationMapIn.at<float>(i,j) <= 1);
            assert(contaminationMapIn.at<float>(i,j) >= 0);

            // we transform the contamination map to a grayscale image to show it
            grayscaleContaminationMap.at<uchar>(i,j) = 255 * ( 1 - contaminationMapIn.at<float>(i,j)); // The darker is the pixel, the more contaminated it is...

        }
    }
    
    // we display the grayscale image
    cv::namedWindow(windowName, cv::WINDOW_AUTOSIZE );
    imshow(windowName, grayscaleContaminationMap);
}


float noiseEstimation(cv::Mat &noisyImageIn, int rMaxIn)
{
    
    assert(noisyImageIn.rows == noisyImageIn.cols); //  the image is squared
    int n = noisyImageIn.rows;
    
    float p_estim = 0;
    float p_estim_1 = 0;
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int nb_black_pixels = 0;
            int nb_pixels = 0;
            for (int i_r = -rMaxIn; i_r <= rMaxIn; i_r++)
            {
                for (int j_r = -rMaxIn; j_r <= rMaxIn; j_r++)
                {
                    if (i + i_r >= 0 && i + i_r < n && j + j_r >= 0 && j + j_r < n) // which means that the pixel (i_r,j_r) is not outside the image
                    {
                        float r = sqrt(float((i_r) ^ 2 + (j_r) ^ 2));
                        if (r <= rMaxIn) // which means that the pixel (i_r,j_r) is in the windows of the filter
                        {
                            nb_pixels++;
                            if (noisyImageIn.at<bool>(i + i_r,j + j_r) == true) // the pixel considered in the windows is black
                            {
                                nb_black_pixels++;
                            }
                        }
                    }
                }
            }
            
            if (nb_black_pixels < (float)nb_pixels / 2)
            {
                p_estim_1 = (float)nb_black_pixels / nb_pixels;
            }
            else
            {
                p_estim_1 = (float)(nb_pixels - nb_black_pixels) / nb_pixels;
            }
            p_estim   = p_estim + p_estim_1;
            
        }
    }
    
    p_estim = (float) p_estim / (n*n);
    return p_estim;
    
}


void computeContaminationAndDistanceMaps(cv::Mat &noisyImageIn, int rMaxIn, float pEstimation, cv::Mat &contaminationMapOut, cv::Mat &distanceMapOut)
{
    assert(noisyImageIn.rows == noisyImageIn.cols);                 //  the image is squared
    assert(contaminationMapOut.rows == contaminationMapOut.cols);   //  the contamination map matrix is squared
    assert(contaminationMapOut.rows == noisyImageIn.rows);          // the image and the contamination map matrix are the same size
    int n = noisyImageIn.rows;
    
    
    // for each pixel in the image...
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int nb_pixels_ligne = 0;
            int nb_pixels = 0;
            int nb_black_pixels = 0;
            std::vector<int> Matrix_nb_pixel_ligne(rMaxIn, 0);
            
            //********************************************************************************************//
            // We count the number of pixels in the window, and the number of pixels per line in the window
            //********************************************************************************************//
            
            
            for (int i_r = -rMaxIn; i_r <= rMaxIn; i_r++)
            {
                nb_pixels_ligne = 0; //Initialization of the number of pixels in a line
                for (int j_r = -rMaxIn; j_r <= rMaxIn; j_r++)
                {
                    if (i + i_r >= 0 && i + i_r < n && j + j_r >= 0 && j + j_r < n) // we make sure that the pixel (i + i_r, j + j_r) is not outside the image
                    {
                        float r = sqrt(float(i_r * i_r + j_r * j_r));
                        if (r <= rMaxIn) // which means that the pixel (i_r,j_r) is in the circle windows of the filter
                        {
                            nb_pixels_ligne++;  // We count the number of pixels per line for implementing the distance
                            nb_pixels++;        // We count the number of pixels in the window for implementing the contamination rate
                            if (noisyImageIn.at<bool>(i + i_r,j + j_r) == true) // the pixel considered in the windows is black
                            {
                                nb_black_pixels++; // We count the number of black pixels in the window for implementing the contamination rate
                            }
                        }
                        
                    }
                }
                
                if (i < (float) n / 2 && i_r > 0) // We just memorize the number of pixels in each line in the half bottom part of the window
                {
                    Matrix_nb_pixel_ligne[i_r - 1] = nb_pixels_ligne;
                }
                
                if (i >= (float) n / 2 && i_r < 0) // We just memorize the number of pixels in each line in the half upper part of the window
                {
                    Matrix_nb_pixel_ligne[i_r + rMaxIn] = nb_pixels_ligne;
                }
                
            }
            
            
            //******************************************//
            // Estimation of the contamination rate :
            //******************************************//
            
            float cont = (float)(nb_pixels * pEstimation - nb_black_pixels) / (nb_pixels * (2 * pEstimation - 1));
            contaminationMapOut.at<float>(i,j) = cont; // contamination for the filter at the pixel (i,j)
            
            if (cont > 1) // Cont_matrix is the contamination rate matrix : each element is < 1 and > 0 in theory.
            {
                contaminationMapOut.at<float>(i,j) = 1;
            }
            else if (cont < 0)
            {
                contaminationMapOut.at<float>(i,j) = 0;
            }
            
            if (cont > (float) 0.5)
            {
                //contaminationMapOut.at<float>(i,j) = 1 - cont;
                contaminationMapOut.at<float>(i,j) = 1 - contaminationMapOut.at<float>(i,j);
            }
            
            
            
            //****************************************************************************************//
            // Estimation of the distance between the center of the filtre and the discountinuity :
            //****************************************************************************************//
            
            // Initialization of sum :
            int sum = 0;
            
            if (i < (float) n / 2)
            {
                int k = rMaxIn;
                while (contaminationMapOut.at<float>(i,j) * nb_pixels > sum && k >= 1)
                {
                    sum = sum + Matrix_nb_pixel_ligne[k - 1];
                    k--;
                }
                float dist = rMaxIn + 0.5 - (rMaxIn - k );
                distanceMapOut.at<float>(i,j) = dist;
            }
            
            if (i >= (float) n / 2)
            {
                int k = 1;
                while (contaminationMapOut.at<float>(i,j) * nb_pixels > sum && k <= rMaxIn)
                {
                    sum = sum + Matrix_nb_pixel_ligne[k - 1];
                    k++;
                }
                float dist = rMaxIn + 0.5 - k + 1;
                distanceMapOut.at<float>(i,j) = dist;
            }
            
            // Computation of the number of pixel per line

            
        }
    }
}


void denoiseFixedCircleMF(cv::Mat &noisyImageIn, int rMaxIn, cv::Mat &denoisedImageOut)
{
    
    // we make sure we have the right inputs
    assert(noisyImageIn.rows == noisyImageIn.cols);                 //  the image is squared
    assert(denoisedImageOut.rows == denoisedImageOut.cols);   //  the contamination map matrix is squared
    assert(denoisedImageOut.rows == noisyImageIn.rows);          // the image and the contamination map matrix are the same size
    int n = noisyImageIn.rows;
    
    
    // for each pixel in the image...
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int nb_pixels = 0;
            int nb_black_pixels = 0;
            
            // We count the number of black pixels in the window
            for (int i_r = -rMaxIn; i_r <= rMaxIn; i_r++)
            {
                for (int j_r = -rMaxIn; j_r <= rMaxIn; j_r++)
                {
                    if (i + i_r >= 0 && i + i_r < n && j + j_r >= 0 && j + j_r < n) // we make sure that the pixel (i + i_r, j + j_r) is not outside the image
                    {
                        float r = sqrt(float(i_r * i_r + j_r * j_r));
                        if (r <= rMaxIn) // which means that the pixel (i_r,j_r) is in the circle windows of the filter
                        {
                            nb_pixels++;        // We count the number of pixels in the window for implementing the contamination rate
                            if (noisyImageIn.at<bool>(i + i_r,j + j_r) == true) // the pixel considered in the windows is black
                            {
                                nb_black_pixels++; // We count the number of black pixels in the window for implementing the contamination rate
                            }
                        }
                        
                    }
                }
            }
            
            // we denoise the image
            if ((float)nb_black_pixels / nb_pixels > 0.5) {
                denoisedImageOut.at<bool>(i,j) = 1; // then pixel is black
            }
            else {
                denoisedImageOut.at<bool>(i,j) = 0; // then pixel is white
            }
            
        }
    }
    
    
}


void denoiseAdaptiveCircleMF(cv::Mat const &noisyImageIn, cv::Mat const &distanceMapIn, int rMaxIn, cv::Mat &denoisedImageOut)
{
    
    // we make sure we have the right inputs
    assert(noisyImageIn.rows == noisyImageIn.cols);                 //  the image is squared
    assert(denoisedImageOut.rows == denoisedImageOut.cols);   //  the contamination map matrix is squared
    assert(denoisedImageOut.rows == noisyImageIn.rows);          // the image and the contamination map matrix are the same size
    int n = noisyImageIn.rows;
    
    
    // for each pixel in the image...
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) {
            
             // we get the estimated distance to a discontinuity
            int part_ent_dist;
            part_ent_dist = floor(distanceMapIn.at<float>(i,j));
            
            // we set the size of the adaptive filter
            float adapt_r = rMaxIn;
            if (part_ent_dist <= (double)rMaxIn / 5)
            {
                adapt_r = (float) rMaxIn / 5; // is the radius at the circle filter if the filter is too close to the discontinuity
            }
            else if (part_ent_dist < (float)rMaxIn / 2)
            {
                adapt_r = part_ent_dist; // is the radius at the circle filter if the filter is relatively close to the discontinuity
            }
            else
            {
                adapt_r = rMaxIn; // is the radius af the circle filter if the filter is far enough  from the discontinuity
            }
            
            // we count the number of black pixels and the number of pixels in the window filter
            int nb_pixels = 0;
            int nb_black_pixels = 0;
            
            for (int i_r = -adapt_r; i_r <= adapt_r; i_r++)
            {
                for (int j_r = -adapt_r; j_r <= adapt_r; j_r++)
                {
                    if (i + i_r >= 0 && i + i_r < n && j + j_r >= 0 && j + j_r < n) // which means that the pixel (i_r,j_r) is not outside the image
                    {
                        
                        float r = sqrt(float(i_r * i_r + j_r * j_r));
                        if (r <= adapt_r) // which means that the pixel (i_r,j_r) is in the windows of the filter
                        {
                            nb_pixels ++; // We count the number of pixels in the window to implement the MF filter
                            
                            if (noisyImageIn.at<bool>(i + i_r,j + j_r) == true) // the pixel considered in the windows is black
                            {
                                nb_black_pixels ++; // We count the number of black pixels in the window for implementing the MF filter
                            }
                        }
                        
                    }
                }
            }
            
            // we set the value of the denoised filter
            if ((float)nb_black_pixels / nb_pixels > 0.5)
            {
                denoisedImageOut.at<bool>(i,j) = true;
            }
            else
            {
                denoisedImageOut.at<bool>(i,j) = false;
            }
            
        }
    }
}




void denoiseAdaptiveEllipticMF(cv::Mat const &noisyImageIn, cv::Mat const &distanceMapIn, float pEstimationIn, int rMaxIn, int numberOfOrientationsIn, cv::Mat &denoisedImageOut)
{
    
    // we make sure we have the right inputs
    assert(noisyImageIn.rows == noisyImageIn.cols);                 //  the image is squared
    assert(denoisedImageOut.rows == denoisedImageOut.cols);   //  the contamination map matrix is squared
    assert(denoisedImageOut.rows == noisyImageIn.rows);          // the image and the contamination map matrix are the same size
    int n = noisyImageIn.rows;
    
    
    // for each pixel in the image...
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) {
            float cont_min = 1; // Initialization of the minimum of the contamination rate (over all the orientations of the ellipses
            int nb_black_pixels_min;
            int nb_pixels_min;
            
            // we get the estimated distance to a discontinuity
            int part_ent_dist;
            part_ent_dist = floor(distanceMapIn.at<float>(i,j));
            
            // we set the small axe length of the adaptive filter
            float adapt_r = rMaxIn;
            if (part_ent_dist <= (double)rMaxIn / 5)
            {
                adapt_r = (float) rMaxIn / 5; // is the radius at the circle filter if the filter is too close to the discontinuity
            }
            else if (part_ent_dist < (float)rMaxIn / 2)
            {
                adapt_r = part_ent_dist; // is the radius at the circle filter if the filter is relatively close to the discontinuity
            }
            else
            {
                adapt_r = rMaxIn; // is the radius af the circle filter if the filter is far enough  from the discontinuity
            }

            
            // for each orientation of the ellipse
            for (int coeff = 0; coeff <= numberOfOrientationsIn - 1; coeff++)
            {
                float theta = (float) coeff / numberOfOrientationsIn;
                theta = theta * CV_PI;
                
                float j_rot = j * cos(theta) - i * sin(theta);
                float i_rot = j * sin(theta) + i * cos(theta);
                
                // we count the number of black pixels and the number of pixels in the window filter
                int nb_pixels = 0;
                int nb_black_pixels = 0;
                for (int i_r = -rMaxIn; i_r <= rMaxIn; i_r++)
                {
                    for (int j_r = -rMaxIn; j_r <= rMaxIn; j_r++)
                    {
                        float j_r_rot = (j + j_r) * cos(theta) - (i + i_r) * sin(theta);
                        float i_r_rot = (j + j_r) * sin(theta) + (i + i_r) * cos(theta);
                        if (i + i_r >= 0 && i + i_r < n && j + j_r >= 0 && j + j_r < n) // which means that the pixel (i_r,j_r) is not outside the image
                        {
                            float r_i = (float) i_rot - i_r_rot ;
                            r_i = (float) r_i / adapt_r;
                            float r_j = (float) j_rot - j_r_rot ;
                            r_j = (float) r_j / rMaxIn;
                            float r = (float) r_i * r_i + r_j * r_j; // equation of an ellipse
                            if (r <= 1) // which means that the pixel (i_r,j_r) is in the windows of the elliptic filter
                            {
                                nb_pixels ++; // We count the number of pixels in the window to implement the MF filter
                                
                                if (noisyImageIn.at<bool>(i + i_r,j + j_r) == true) // the pixel considered in the windows is black
                                {
                                    nb_black_pixels ++; // We count the number of black pixels in the window for implementing the MF filter
                                }
                            }
                            
                        }
                    }
                }
                
                // we compute the contamination rate of the curent position of the ellipse
                float cont_e_num = (nb_pixels * pEstimationIn - nb_black_pixels);
                float cont_e_den = nb_pixels * (2 * pEstimationIn - 1);
                float cont_e = (double)cont_e_num / cont_e_den; // contamination rate in the ellipse;
                
                // we make sure the contamination rate has valu that makes sense
                if (cont_e > 1) // cont is the contamination rate matrix : each element is < 1 and > 0 in theory.
                {
                    cont_e = 1;
                }
                else if (cont_e < 0)
                {
                    cont_e = 0;
                }
                if (cont_e > 0.5)
                {
                    cont_e = 1 - cont_e;
                }
                
                // we check out if the contamination rate is smaller than the minimum contamination rate (which is stored)
                if (cont_e <= cont_min)
                {
                    cont_min = cont_e;                      // if so, we update the minimum contamination rate
                    nb_black_pixels_min = nb_black_pixels;  // we store the number of black pixels for that configuration
                    nb_pixels_min = nb_pixels;              // we store the number of pixels for that configuration
                }
            }
            
            
            // we set the value of the denoised filter
            if ((float)nb_black_pixels_min / nb_pixels_min > 0.5)
            {
                denoisedImageOut.at<bool>(i,j) = true;
            }
            else
            {
                denoisedImageOut.at<bool>(i,j) = false;
            }

            
        }
    }
}


float computeErrorRate(cv::Mat const &image1In, cv::Mat const &image2In)
{
    // we make sure we have the right inputs
    assert(image1In.rows == image2In.cols);   //  the image 1 is squared
    assert(image1In.rows == image2In.rows);   //  image1 and image2 have to be same size
    assert(image1In.rows == image2In.cols);   //  image1 and image2 have to be same size
    int n = image1In.rows;
    
    
    // we compute the error rate
    int errorCounter = 0;
    for (int i = 0; i < n; i++) // for each pixel in the image...
    {
        for (int j = 0; j < n; j++) {
            if (image1In.at<bool>(i,j) != image2In.at<bool>(i,j))
            {
                errorCounter ++; // we increment the error counter
            }
        }
    }
    float errorRate = errorCounter / (n*n);
    
    // we return the error rate
    return errorRate;
    
}






