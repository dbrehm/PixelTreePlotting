#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"

// This code is adapted from previous code. It is NOT optimized as I was adding things ad hoc as I went.
// It also uses the PixelTrees data format so you should go here: https://github.com/CMSTrackerDPG/SiPixelTools-PixelTrees/blob/master/plugins/PixelTree.cc for what is accessible in the TTree
// The overall structure is I have made two functions to fill histograms. Some computing is done for each.
// Then I run over the TProfile of the Y Pulls to get the RMS and error of RMS and fill a 1d Histogram of that.
// Then at the plotting step things are fit, including the Y Pull RMS vs. ClCharge/ClSize fit.


// NOTE!!!!!
// MUST MANUALLY CHANGE NAMING SCHEMES, LAYER OR DISK BEING SELECTED, ERROR CORRECTION PROFILE, ETC....

void fill_pulls(std::string iFile, TH1F *h_xpullrms_clX, TH1F *h_ypullrms_clY, TH1F *h_x, TH1F *h_y, TH1F *h_xe, TH1F *h_ye, TH1F *h_ClCheck_x, TH1F *h_ClCheck_y, TH1F *h_xrms_clX, TH1F *h_yrms_clY, TProfile *h_x_clX, TProfile *h_y_clY, TProfile *h_xe_clX, TProfile *h_ye_clY, TH2F *h_xpullrms_clX_qbin, TH2F *h_ypullrms_clY_qbin, TH1F *h_pullx = nullptr, TH1F *h_pully = nullptr, TProfile *h_pullx_clX = nullptr, TProfile *h_pully_clY = nullptr, bool on_edge = false, bool bad_pix = false, bool make_3d = false, bool pass2 = false)
{

    TFile *f_data = TFile::Open(iFile.c_str());
    TTree *t1 = (TTree *)f_data->Get("pixelTree");

    //read event data
    const int SIZE = 20000;
    Long64_t size  =  t1->GetEntries();
    int count = 0;
    Float_t TkEta[10000];
    Float_t TkPt[10000];
    Float_t ClRhLx[SIZE], ClRhLy[SIZE], ClRhLxE[SIZE], ClRhLyE[SIZE], ClCharge[SIZE], ClSimTrEta[SIZE][10], ClSimHitLx[SIZE][10], ClSimHitLy[SIZE][10];
    Int_t ClRhIsOnEdge[SIZE], ClN, ClSimHitN[SIZE], ClType[SIZE], ClRhHasBadPix[SIZE], ClSize[SIZE], ClSizeX[SIZE], ClSizeY[SIZE], ClLayer[SIZE], ClDisk[SIZE], ClRhqBin[SIZE];

    t1->SetBranchAddress("TkPt", &TkPt);
    t1->SetBranchAddress("TkEta", &TkEta);
    t1->SetBranchAddress("ClN", &ClN);
    t1->SetBranchAddress("ClSimHitN", &ClSimHitN);
    t1->SetBranchAddress("ClType", &ClType);
    t1->SetBranchAddress("ClRhHasBadPixels", &ClRhHasBadPix);
    t1->SetBranchAddress("ClRhIsOnEdge", &ClRhIsOnEdge);
    t1->SetBranchAddress("ClRhLx", &ClRhLx);
    t1->SetBranchAddress("ClRhLy", &ClRhLy);
    t1->SetBranchAddress("ClRhLxE", &ClRhLxE);
    t1->SetBranchAddress("ClRhLyE", &ClRhLyE);
    t1->SetBranchAddress("ClSimTrEta", &ClSimTrEta);
    t1->SetBranchAddress("ClSimHitLx", &ClSimHitLx);
    t1->SetBranchAddress("ClSimHitLy", &ClSimHitLy);
    t1->SetBranchAddress("ClSize", &ClSize);
    t1->SetBranchAddress("ClSizeX", &ClSizeX);
    t1->SetBranchAddress("ClSizeY", &ClSizeY);
    t1->SetBranchAddress("ClLayer", &ClLayer);
    t1->SetBranchAddress("ClDisk", &ClDisk);
    t1->SetBranchAddress("ClRhqBin", &ClRhqBin);
    t1->SetBranchAddress("ClCharge", &ClCharge);

    for(int i =0; i< size; i++){
        t1->GetEntry(i);
        if(ClN > SIZE) printf("WARNING: Number of clusters (%i) greater than array size (%i), memory issues likely!!! \n \n \n", ClN, SIZE);
        for(int j=0; j<ClN; j++){
            // HERE YOU CHANGE CLayer TO THE ONE YOU WANT
            if(ClType[j]==1  && ClLayer[j]== 1)// && ClRhqBin[j]==0) // on a track
            {    // std::cout << "clust layer " << ClLayer[j]  << std::endl
                //std::cout << "tk eta " << fabs(TkEta[0]) << " sim eta " << fabs(ClSimTrEta[j][0]) << std::endl;
                float dx(9999), dy(9999), pullx(9999), pully(9999);
                int iKx(0), iKy(0);
                for(int k=0; k<ClSimHitN[j]; k++){
                    if(fabs(ClSimHitLx[j][k] - ClRhLx[j]) < fabs(dx)){ 
                        dx = ClRhLx[j] - ClSimHitLx[j][k];
                        pullx =  dx/ClRhLxE[j];
                        iKx = k;
                    }
                    if(fabs(ClSimHitLy[j][k] - ClRhLy[j]) < fabs(dy)){ 
                        dy =  (ClRhLy[j]  - ClSimHitLy[j][k]);
                        // HERE YOU CHANGE ERROR CORRECTION SCHEME
                        // pully = dy /ClRhLyE[j];
                        pully = dy / (ClRhLyE[j] * (0.0003 * (pow(ClCharge[j] / ClSizeY[j], 2) - 0.0045 * pow(ClCharge[j] / ClSizeY[j], 1)) + 1.656));

                        iKy = k;
                    }
                }
                bool fill = (dx<9999) && (dy<9999);
                if(on_edge) fill = fill && ClRhIsOnEdge[j];
                if(bad_pix) fill = fill && ClRhHasBadPix[j];
                if(fill)// && (abs(pully) <= 1.0))
                {
                    //mult by 10000 to convert to microns
                    float to_microns = 1e4;
                    if(make_3d == false && pass2 == false){
                        h_ClCheck_x->Fill(ClSizeX[j]);
                        h_ClCheck_y->Fill(ClSizeY[j]);

                        h_x->Fill(to_microns * dx);
                        h_y->Fill(to_microns * dy);

                        h_xe->Fill(to_microns * ClRhLxE[j]);
                        h_ye->Fill(to_microns * ClRhLyE[j]);

                        // std::cout << " check = " << ClCharge[j] / ClSize[j] << std::endl;
                        // if (ClCharge[j] / ClSize[j] > 22.0){
                        //     count++;
                        //     std::cout << "clust Charge is " << ClCharge[j] << " clust size " << ClSize[j] << " clust size x is " << ClSizeX[j] << "clust size y is " << ClSizeY[j] << std::endl;
                        //     std::cout << "The count is: " << count << std::endl;
                        // }
                        // std::cout << " clust size = " << ClSize[j] << " clust size x * clust size y = " << ClSizeX[j] * ClSizeY[j] << std::endl;
                        h_xe_clX->Fill(ClCharge[j] / ClSizeX[j], to_microns * ClRhLxE[j], 1);
                        h_ye_clY->Fill(ClCharge[j] / ClSizeY[j], to_microns * ClRhLyE[j], 1);

                        h_x_clX->Fill(ClSizeX[j], to_microns * dx, 1);
                        h_y_clY->Fill(ClSizeY[j], to_microns * dy, 1);
                    }
                    if (pass2 == true){
                        // int xbin = h_x_clX->FindBin(ClSizeX[j]);
                        // int ybin = h_y_clY->FindBin(ClSizeY[j]);
                        // h_xrms_clX->Fill(ClSizeX[j], h_x_clX->GetBinError(xbin), 1);
                        // h_yrms_clY->Fill(ClSizeY[j], h_y_clY->GetBinError(ybin), 1);

                        // int xPbin = h_pullx_clX->FindBin(ClSizeX[j]);
                        // int yPbin = h_pully_clY->FindBin(ClSizeY[j]);
                        // h_xpullrms_clX->Fill(ClSizeX[j], h_pullx_clX->GetBinError(xPbin), 1);
                        // h_ypullrms_clY->Fill(ClSizeY[j], h_pully_clY->GetBinError(yPbin), 1);
                    }else{
                        if (h_pullx != nullptr)
                        {

                            h_xpullrms_clX_qbin->Fill(ClSizeX[j], pullx);
                            // int pull_bin = h_pullx_clX->FindBin(ClSizeX[j],pullx);
                            // h_xpullrms_clX_qbin->SetBinContent(bin, ClRhqBin[j]);

                            h_pullx->Fill(pullx);
                            h_pullx_clX->Fill(ClCharge[j] / ClSizeX[j], pullx, 1);
                        }
                        if (h_pully != nullptr)
                        {

                            int bin = h_ypullrms_clY_qbin->Fill(ClSizeY[j], pully);
                            // int pull_bin = h_pully_clY->FindBin(ClSizeY[j],pully);
                            // h_ypullrms_clY_qbin->SetBinContent(bin, ClRhqBin[j]);
                            // if (abs(pully) > 1.0)
                            // {
                            //     pully = (1 / (0.1*(log(ClCharge[j] / ClSizeY[j])))) * (dy / ClRhLyE[j]);
                            // }
                            h_pully->Fill(pully);
                            h_pully_clY->Fill(ClCharge[j] / ClSizeY[j], pully, 1);
                        }
                    }

                    
                }
            }
        }
    }
    printf("printing means and std devs for %s (N = %.0f) \n", iFile.c_str(), h_x->Integral());
    if(on_edge) printf("Edge Clusters: \n");
    if(bad_pix) printf("BadPix Clusters: \n");
    printf("X: mean %.3e std_dev %.3e \n", h_x->GetMean(1), h_x->GetStdDev(1));
    printf("Y: mean %.3e std_dev %.3e \n", h_y->GetMean(1), h_y->GetStdDev(1));
    if(h_pullx != nullptr) printf("Pull X: mean %.3e std_dev %.3e \n", h_pullx->GetMean(1), h_pullx->GetStdDev(1));
    if(h_pully != nullptr) printf("Pull Y: mean %.3e std_dev %.3e \n", h_pully->GetMean(1), h_pully->GetStdDev(1));

    //h_x->Scale(1/h_x->Integral());
    //h_y->Scale(1/h_y->Integral());

    //if(h_pullx != nullptr) h_pullx->Scale(1/h_pullx->Integral());
    //if(h_pully != nullptr) h_pully->Scale(1/h_pully->Integral());

    return;
}

void fill_pulls24(std::string iFile, TH1F *h_xpullrms_clX_24, TH1F *h_ypullrms_clY_24, TH1F *h_x_24, TH1F *h_y_24, TH1F *h_xe_24, TH1F *h_ye_24, TH1F *h_ClCheck_x_24, TH1F *h_ClCheck_y_24, TH1F *h_xrms_clX_24, TH1F *h_yrms_clY_24, TProfile *h_x_clX_24, TProfile *h_y_clY_24, TProfile *h_xe_clX_24, TProfile *h_ye_clY_24, TH2F *h_xpullrms_clX_qbin_24, TH2F *h_ypullrms_clY_qbin_24, TH1F *h_pullx_24 = nullptr, TH1F *h_pully_24 = nullptr, TProfile *h_pullx_clX_24 = nullptr, TProfile *h_pully_clY_24 = nullptr, bool on_edge = false, bool bad_pix = false, bool make_3d = false, bool pass2 = false)
{
    // TFile *f_data = TFile::Open(iFile.c_str());
    TChain t1("pixelTree");
    t1.Add("PixelTree_wDetAngles_1.root");
    t1.Add("PixelTree_wDetAngles_2.root");
    t1.Add("PixelTree_wDetAngles_3.root");
    t1.Add("PixelTree_wDetAngles_4.root");
    t1.Add("PixelTree_wDetAngles_5.root");

    //read event data
    const int SIZE = 20000;
    Long64_t size = t1.GetEntries();
    int count = 0;
    Float_t TkEta[10000];
    Float_t TkPt[10000];
    Float_t ClRhLx[SIZE], ClRhLy[SIZE], ClRhLxE[SIZE], ClRhLyE[SIZE], ClCharge[SIZE], ClSimTrEta[SIZE][10], ClSimHitLx[SIZE][10], ClSimHitLy[SIZE][10];
    Int_t ClRhIsOnEdge[SIZE], ClN, ClSimHitN[SIZE], ClType[SIZE], ClRhHasBadPix[SIZE], ClSize[SIZE], ClSizeX[SIZE], ClSizeY[SIZE], ClLayer[SIZE], ClDisk[SIZE], ClRhqBin[SIZE];

    t1.SetBranchAddress("TkPt", &TkPt);
    t1.SetBranchAddress("TkEta", &TkEta);
    t1.SetBranchAddress("ClN", &ClN);
    t1.SetBranchAddress("ClSimHitN", &ClSimHitN);
    t1.SetBranchAddress("ClType", &ClType);
    t1.SetBranchAddress("ClRhHasBadPixels", &ClRhHasBadPix);
    t1.SetBranchAddress("ClRhIsOnEdge", &ClRhIsOnEdge);
    t1.SetBranchAddress("ClRhLx", &ClRhLx);
    t1.SetBranchAddress("ClRhLy", &ClRhLy);
    t1.SetBranchAddress("ClRhLxE", &ClRhLxE);
    t1.SetBranchAddress("ClRhLyE", &ClRhLyE);
    t1.SetBranchAddress("ClSimTrEta", &ClSimTrEta);
    t1.SetBranchAddress("ClSimHitLx", &ClSimHitLx);
    t1.SetBranchAddress("ClSimHitLy", &ClSimHitLy);
    t1.SetBranchAddress("ClSize", &ClSize);
    t1.SetBranchAddress("ClSizeX", &ClSizeX);
    t1.SetBranchAddress("ClSizeY", &ClSizeY);
    t1.SetBranchAddress("ClLayer", &ClLayer);
    t1.SetBranchAddress("ClDisk", &ClDisk);
    t1.SetBranchAddress("ClRhqBin", &ClRhqBin);
    t1.SetBranchAddress("ClCharge", &ClCharge);

    for (int i = 0; i < size; i++)
    {
        t1.GetEntry(i);
        if (ClN > SIZE)
            printf("WARNING: Number of clusters (%i) greater than array size (%i), memory issues likely!!! \n \n \n", ClN, SIZE);
        for (int j = 0; j < ClN; j++)
        {
            if (ClType[j] == 1 && ClLayer[j] == 1) // && ClRhqBin[j]==0)
            { // on a track
                // std::cout << "clust layer " << ClLayer[j]  << std::endl
                //std::cout << "tk eta " << fabs(TkEta[0]) << " sim eta " << fabs(ClSimTrEta[j][0]) << std::endl;
                float dx(9999), dy(9999), pullx(9999), pully(9999);
                int iKx(0), iKy(0);
                for (int k = 0; k < ClSimHitN[j]; k++)
                {
                    if (fabs(ClSimHitLx[j][k] - ClRhLx[j]) < fabs(dx))
                    {
                        dx = ClRhLx[j] - ClSimHitLx[j][k];
                        pullx = dx / ClRhLxE[j];
                        iKx = k;
                    }
                    if (fabs(ClSimHitLy[j][k] - ClRhLy[j]) < fabs(dy))
                    {
                        dy = (ClRhLy[j] - ClSimHitLy[j][k]);
                        // pully = dy / ClRhLyE[j];
                        pully = dy / (ClRhLyE[j]*(0.00028*(pow(ClCharge[j] / ClSizeY[j], 2) - 0.025*pow(ClCharge[j] / ClSizeY[j], 1)) + 2.15));

                        iKy = k;
                    }
                }
                bool fill = (dx < 9999) && (dy < 9999);
                if (on_edge)
                    fill = fill && ClRhIsOnEdge[j];
                if (bad_pix)
                    fill = fill && ClRhHasBadPix[j];
                if (fill)// && (abs(pully) <= 1.0))
                {
                    //mult by 10000 to convert to microns
                    float to_microns = 1e4;
                    if (make_3d == false && pass2 == false)
                    {
                        h_ClCheck_x_24->Fill(ClSizeX[j]);
                        h_ClCheck_y_24->Fill(ClSizeY[j]);

                        h_x_24->Fill(to_microns * dx);
                        h_y_24->Fill(to_microns * dy);

                        h_xe_24->Fill(to_microns * ClRhLxE[j]);
                        h_ye_24->Fill(to_microns * ClRhLyE[j]);

                        // std::cout << " check = " << ClCharge[j] / ClSize[j] << std::endl;
                        // if (ClCharge[j] / ClSize[j] > 22.0){
                        //     count++;
                        //     std::cout << "clust Charge is " << ClCharge[j] << " clust size " << ClSize[j] << " clust size x is " << ClSizeX[j] << "clust size y is " << ClSizeY[j] << std::endl;
                        //     std::cout << "The count is: " << count << std::endl;
                        // }
                        // std::cout << " clust size = " << ClSize[j] << " clust size x * clust size y = " << ClSizeX[j] * ClSizeY[j] << std::endl;
                        h_xe_clX_24->Fill(ClCharge[j] / ClSizeX[j], to_microns * ClRhLxE[j], 1);
                        h_ye_clY_24->Fill(ClCharge[j] / ClSizeY[j], to_microns * ClRhLyE[j], 1);

                        h_x_clX_24->Fill(ClSizeX[j], to_microns * dx, 1);
                        h_y_clY_24->Fill(ClSizeY[j], to_microns * dy, 1);
                    }
                    if (pass2 == true)
                    {
                        // int xbin = h_x_clX_24->FindBin(ClSizeX[j]);
                        // int ybin = h_y_clY_24->FindBin(ClSizeY[j]);
                        // h_xrms_clX_24->Fill(ClSizeX[j], h_x_clX_24->GetBinError(xbin), 1);
                        // h_yrms_clY_24->Fill(ClSizeY[j], h_y_clY_24->GetBinError(ybin), 1);

                        // int xPbin = h_pullx_clX_24->FindBin(ClSizeX[j]);
                        // int yPbin = h_pully_clY_24->FindBin(ClSizeY[j]);
                        // h_xpullrms_clX_24->Fill(ClSizeX[j], h_pullx_clX_24->GetBinError(xPbin), 1);
                        // h_ypullrms_clY_24->Fill(ClSizeY[j], h_pully_clY_24->GetBinError(yPbin), 1);
                    }else{
                        if (h_pullx_24 != nullptr)
                        {

                            h_xpullrms_clX_qbin_24->Fill(ClSizeX[j], pullx);
                            // int pull_bin = h_pullx_clX->FindBin(ClSizeX[j],pullx);
                            // h_xpullrms_clX_qbin->SetBinContent(bin, ClRhqBin[j]);

                            h_pullx_24->Fill(pullx);
                            h_pullx_clX_24->Fill(ClCharge[j] / ClSizeX[j], pullx, 1);
                        }
                        if (h_pully_24 != nullptr)
                        {

                            int bin = h_ypullrms_clY_qbin_24->Fill(ClSizeY[j], pully);
                            // int pull_bin = h_pully_clY->FindBin(ClSizeY[j],pully);
                            // h_ypullrms_clY_qbin->SetBinContent(bin, ClRhqBin[j]);

                            // if(abs(pully) > 1.0){
                            //     pully = (1 / (0.1*(log(ClCharge[j] / ClSizeY[j])))) * (dy / ClRhLyE[j]);
                            // }
                            h_pully_24->Fill(pully);
                            h_pully_clY_24->Fill(ClCharge[j] / ClSizeY[j], pully, 1);
                        }
                    }
                }
            }
        }
    }
    printf("printing means and std devs for %s (N = %.0f) \n", iFile.c_str(), h_x_24->Integral());
    if (on_edge)
        printf("Edge Clusters: \n");
    if (bad_pix)
        printf("BadPix Clusters: \n");
    printf("X: mean %.3e std_dev %.3e \n", h_x_24->GetMean(1), h_x_24->GetStdDev(1));
    printf("Y: mean %.3e std_dev %.3e \n", h_y_24->GetMean(1), h_y_24->GetStdDev(1));
    if (h_pullx_24 != nullptr)
        printf("Pull X: mean %.3e std_dev %.3e \n", h_pullx_24->GetMean(1), h_pullx_24->GetStdDev(1));
    if (h_pully_24 != nullptr)
        printf("Pull Y: mean %.3e std_dev %.3e \n", h_pully_24->GetMean(1), h_pully_24->GetStdDev(1));

    //h_x->Scale(1/h_x->Integral());
    //h_y->Scale(1/h_y->Integral());

    //if(h_pullx != nullptr) h_pullx->Scale(1/h_pullx->Integral());
    //if(h_pully != nullptr) h_pully->Scale(1/h_pully->Integral());

    return;
}
void draw_pulls(){
    float range = 300.;
    float yrange = 80;
    float range2 = 10;
    int nBins = 100;
    int nYBins = 20;
    int bins2 = 10;

    gROOT->SetBatch(1);
    gStyle->SetOptFit(1);
    gStyle->SetPalette(55);

    TH1F *h_ClCheck_x = new TH1F("h_ClCheck_x", "X Cluster Check 2021_layer1_fitChargePerSizeSkipPol; Cluster Size X (Pixels)", nYBins, 0.5, 20.5);
    TH1F *h_ClCheck_y = new TH1F("h_ClCheck_y", "Y Cluster Check 2021_layer1_fitChargePerSizeSkipPol; Cluster Size Y (Pixels)", nYBins, 0.5, 20.5);

    TH1F *h_x = new TH1F("h_x", "X position residual 2021_layer1_fitChargePerSizeSkipPol; #Deltax (#mum)" , nBins, -range, range);
    TH1F *h_y = new TH1F("h_y", "Y position residual 2021_layer1_fitChargePerSizeSkipPol; #Deltay (#mum)" , nBins, -range, range);

    TH1F *h_xe = new TH1F("h_xe", "X position error 2021_layer1_fitChargePerSizeSkipPol; x Error (#mum)", nBins, 0, 200);
    TH1F *h_ye = new TH1F("h_ye", "Y position error 2021_layer1_fitChargePerSizeSkipPol; y Error (#mum)", nBins, 0, 200);

    TH1F *h_pullx = new TH1F("h_pullx", "X Pull 2021_layer1_fitChargePerSizeSkipPol; X Pull" , nBins,-5, 5);
    TH1F *h_pully = new TH1F("h_pully", "Y Pull 2021_layer1_fitChargePerSizeSkipPol; Y Pull" , nBins,-5, 5);

    TH2F *h_xpullrms_clX_qbin = new TH2F("h_xpullrms_clX_qbin", "qbin0 as X pull RMS v. Cluster Size X 2021_layer1_fitChargePerSizeSkipPol; Cluster Size X (Pixel); x Pull RMS (#mum)" , bins2, 1, 10, bins2, 0, range2);
    TH2F *h_ypullrms_clY_qbin = new TH2F("h_ypullrms_clY_qbin", "qbin0 as Y pull RMS v. Cluster Size Y 2021_layer1_fitChargePerSizeSkipPol; Cluster Size Y (Pixel); y Pull RMS (#mum)", bins2, 1, 10, bins2, 0, range2);

    TProfile *h_x_clX = new TProfile("h_x_clX", "X position residual v. Cluster Size 2021_layer1_fitChargePerSizeSkipPol; Cluster Size (Pixels); #Deltax (#mum);", bins2, 0.5, 10.5, 0, range, "S");
    TProfile *h_y_clY = new TProfile("h_y_clY", "Y position residual v. Cluster Size 2021_layer1_fitChargePerSizeSkipPol; Cluster Size (Pixels); #Deltay (#mum);", bins2, 0.5, 10.5, 0, range, "S");

    TProfile *h_pullx_clX = new TProfile("h_pullx_clX", "X Pull v. Cluster Size 2021_layer1_fitChargePerSizeSkipPol; Cluster Size X(Pixels); X Pull", bins2, 0.5, 30.5,-5, 5, "S");
    TProfile *h_pully_clY = new TProfile("h_pully_clY", "Y Pull v. Cluster Size 2021_layer1_fitChargePerSizeSkipPol; Cluster Size Y(Pixels); Y Pull", bins2, 0.5, 30.5,-5, 5, "S");

    int vBinsNum = 8;
    double vbins[9] = {0.5,3.5,6.5,9.5,12.5,15.5,18.5,21.5,24.5};
    
    TH1F *h_xrms_clX = new TH1F("h_xrms_clX", "X position RMS v. Cluster Size 2021_layer1_fitChargePerSizeSkipPol; Cluster Size X(Pixels); X Res RMS (#mum)", vBinsNum, vbins);
    TH1F *h_yrms_clY = new TH1F("h_yrms_clY", "Y position RMS v. Cluster Size 2021_layer1_fitChargePerSizeSkipPol; Cluster Size Y(Pixels); Y Res RMS (#mum)", vBinsNum, vbins);

    TH1F *h_xpullrms_clX = new TH1F("h_xpullrms_clX", "X Pull RMS v. Cluster Size 2021_layer1_fitChargePerSizeSkipPol;Charge /  Cluster Size X(ke/Pixels); X Pull RMS", vBinsNum, vbins);
    TH1F *h_ypullrms_clY = new TH1F("h_ypullrms_clY", "Y Pull RMS v. Cluster Size 2021_layer1_fitChargePerSizeSkipPol; Charge / Cluster Size Y(ke/Pixels); Y Pull RMS", vBinsNum, vbins);

    TProfile *h_xe_clX = new TProfile("h_xe_clX", "X position error v. Cluster Size layer1; Charge / Cluster Size X (ke/Pixels); x Error (#mum)", 8, 1, yrange, 0, range, "S");
    TProfile *h_ye_clY = new TProfile("h_ye_clY", "Y position error v. Cluster Size layer1; Charge / Cluster Size Y (ke/Pixels); y Error (#mum)", 8, 1, yrange, 0, range, "S");



    h_x->SetMarkerColor(kBlack);
    h_x->SetMarkerStyle(20);
    h_x->SetLineColor(kBlack);

    h_y->SetMarkerColor(kBlack);
    h_y->SetMarkerStyle(20);
    h_y->SetLineColor(kBlack);

    h_pullx->SetMarkerColor(kBlack);
    h_pullx->SetMarkerStyle(20);
    h_pullx->SetLineColor(kBlack);

    h_pully->SetMarkerColor(kBlack);
    h_pully->SetMarkerStyle(20);
    h_pully->SetLineColor(kBlack);

    TH1F *h_ClCheck_x_24 = new TH1F("h_ClCheck_x_24", "X Cluster Check 2024_layer1_fitChargePerSizeSkipPol; Cluster Size X (Pixels)", nYBins, 0.5, 20.5);
    TH1F *h_ClCheck_y_24 = new TH1F("h_ClCheck_y_24", "Y Cluster Check 2024_layer1_fitChargePerSizeSkipPol; Cluster Size Y (Pixels)", nYBins, 0.5, 20.5);

    TH1F *h_x_24 = new TH1F("h_x_24", "X position residual 2024_layer1_fitChargePerSizeSkipPol; #Deltax (#mum)", nBins, -range, range);
    TH1F *h_y_24 = new TH1F("h_y_24", "Y position residual 2024_layer1_fitChargePerSizeSkipPol; #Deltay (#mum)", nBins, -range, range);

    TH1F *h_xe_24 = new TH1F("h_xe_24", "X position error 2024_layer1_fitChargePerSizeSkipPol; x Error (#mum)", nBins, 0, 200);
    TH1F *h_ye_24 = new TH1F("h_ye_24", "Y position error 2024_layer1_fitChargePerSizeSkipPol; y Error (#mum)", nBins, 0, 200);

    TH1F *h_pullx_24 = new TH1F("h_pullx_24", "X Pull 2024_layer1_fitChargePerSizeSkipPol; X Pull", nBins,-5, 5);
    TH1F *h_pully_24 = new TH1F("h_pully_24", "Y Pull 2024_layer1_fitChargePerSizeSkipPol; Y Pull", nBins,-5, 5);

    TH2F *h_xpullrms_clX_qbin_24 = new TH2F("h_xpullrms_clX_qbin", "qbin0 as X pull RMS v. Cluster Size X 2021_layer1_fitChargePerSizeSkipPol; Cluster Size X (Pixel); x Pull RMS (#mum)", bins2, 1, 10, bins2, 0, range2);
    TH2F *h_ypullrms_clY_qbin_24 = new TH2F("h_ypullrms_clY_qbin", "qbin0 as Y pull RMS v. Cluster Size Y 2021_layer1_fitChargePerSizeSkipPol; Cluster Size Y (Pixel); y Pull RMS (#mum)", bins2, 1, 10, bins2, 0, range2);

    TProfile *h_x_clX_24 = new TProfile("h_x_clX_24", "X position residual v. Cluster Size 2024_layer1_fitChargePerSizeSkipPol; Cluster Size (Pixels); #Deltax (#mum);", 40, 1, yrange, 0, range, "S");
    TProfile *h_y_clY_24 = new TProfile("h_y_clY_24", "Y position residual v. Cluster Size 2024_layer1_fitChargePerSizeSkipPol; Cluster Size (Pixels); #Deltay (#mum);", 40, 1, yrange, 0, range, "S");

    TH1F *h_xrms_clX_24 = new TH1F("h_xrms_clX_24", "X position RMS v. Cluster Size 2024_layer1_fitChargePerSizeSkipPol; Cluster Size X(Pixels); X Res RMS (#mum)", vBinsNum, vbins);
    TH1F *h_yrms_clY_24 = new TH1F("h_yrms_clY_24", "Y position RMS v. Cluster Size 2024_layer1_fitChargePerSizeSkipPol; Cluster Size Y(Pixels); Y Res RMS (#mum)", vBinsNum, vbins);

    TH1F *h_xpullrms_clX_24 = new TH1F("h_xpullrms_clX_24", "X Pull RMS v. Cluster Size 2024_layer1_fitChargePerSizeSkipPol; Charge /Cluster Size X(ke/Pixels); X Pull RMS", vBinsNum, vbins);
    TH1F *h_ypullrms_clY_24 = new TH1F("h_ypullrms_clY_24", "Y Pull RMS v. Cluster Size 2024_layer1_fitChargePerSizeSkipPol; Charge /Cluster Size Y(ke/Pixels); Y Pull RMS", vBinsNum, vbins);

    TProfile *h_xe_clX_24 = new TProfile("h_xe_clX_24", "X position error v. Cluster Size 2024_layer1_fitChargePerSizeSkipPol; Charge / Cluster Size X (ke/Pixels); x Error (#mum)", 8, 1, yrange, 0, range, "S");
    TProfile *h_ye_clY_24 = new TProfile("h_ye_clY_24", "Y position error v. Cluster Size 2024_layer1_fitChargePerSizeSkipPol; Charge / Cluster Size Y (ke/Pixels); y Error (#mum)", 8, 1, yrange, 0, range, "S");

    TProfile *h_pullx_clX_24 = new TProfile("h_pullx_clX_24", "X Pull v. Cluster Size layer1; Cluster Size X(Pixels); X Pull", bins2, 0.5, 30.5,-5, 5, "S");
    TProfile *h_pully_clY_24 = new TProfile("h_pully_clY_24", "Y Pull v. Cluster Size layer1; Cluster Size Y(Pixels); Y Pull", bins2, 0.5, 30.5,-5, 5, "S");

    h_x_24->SetMarkerColor(kBlack);
    h_x_24->SetMarkerStyle(20);
    h_x_24->SetLineColor(kBlack);

    h_y_24->SetMarkerColor(kBlack);
    h_y_24->SetMarkerStyle(20);
    h_y_24->SetLineColor(kBlack);

    h_pullx_24->SetMarkerColor(kBlack);
    h_pullx_24->SetMarkerStyle(20);
    h_pullx_24->SetLineColor(kBlack);

    h_pully_24->SetMarkerColor(kBlack);
    h_pully_24->SetMarkerStyle(20);
    h_pully_24->SetLineColor(kBlack);

    string f_name("PixelTree_wDetAngles_2021.root");
    string f_name24("PixelTree_wDetAngles.root");
    //All clusters

    bool use_edge = false;
    bool use_badpix = false;
    bool make_3d = false;
    bool pass2 = false;

    // Fill pulls based on functions above.
    fill_pulls(f_name, h_xpullrms_clX, h_ypullrms_clY, h_x, h_y, h_xe, h_ye, h_ClCheck_x, h_ClCheck_y, h_xrms_clX, h_yrms_clY, h_x_clX, h_y_clY, h_xe_clX, h_ye_clY, h_xpullrms_clX_qbin, h_ypullrms_clY_qbin, h_pullx, h_pully, h_pullx_clX, h_pully_clY, use_edge, use_badpix, make_3d, pass2);

    fill_pulls24(f_name24, h_xpullrms_clX_24, h_ypullrms_clY_24, h_x_24, h_y_24, h_xe_24, h_ye_24, h_ClCheck_x_24, h_ClCheck_y_24, h_xrms_clX_24, h_yrms_clY_24, h_x_clX_24, h_y_clY_24, h_xe_clX_24, h_ye_clY_24, h_xpullrms_clX_qbin_24, h_ypullrms_clY_qbin_24, h_pullx_24, h_pully_24, h_pullx_clX_24, h_pully_clY_24, use_edge, use_badpix, make_3d, pass2);

    // After making the TProfile for the Y Pulls to get RMS we then fill the 1d Histogram with the RMS and error of RMS from the TProfile

    for(int i = 2; i < 8; i++){
        
        h_xrms_clX->SetBinContent(i+1, h_x_clX->GetBinError(i+1));
        h_xrms_clX->SetBinError(i+1, h_x_clX->GetBinError(i+1)/TMath::Sqrt(h_x_clX->GetBinEntries(i+1)));

        h_yrms_clY->SetBinContent(i+1, h_y_clY->GetBinError(i+1));
        h_yrms_clY->SetBinError(i+1, h_y_clY->GetBinError(i+1)/TMath::Sqrt(h_y_clY->GetBinEntries(i+1)));

        h_xpullrms_clX->SetBinContent(i+1, h_pullx_clX->GetBinError(i+1));
        h_xpullrms_clX->SetBinError(i+1, h_pullx_clX->GetBinError(i+1)/TMath::Sqrt(h_pullx_clX->GetBinEntries(i+1)));

        h_ypullrms_clY->SetBinContent(i+1, h_pully_clY->GetBinError(i+1));
        h_ypullrms_clY->SetBinError(i+1, h_pully_clY->GetBinError(i+1)/TMath::Sqrt(h_pully_clY->GetBinEntries(i+1)));



        h_xrms_clX_24->SetBinContent(i+1, h_x_clX_24->GetBinError(i+1));
        h_xrms_clX_24->SetBinError(i+1, h_x_clX_24->GetBinError(i+1)/TMath::Sqrt(h_x_clX_24->GetBinEntries(i+1)));

        h_yrms_clY_24->SetBinContent(i+1, h_y_clY_24->GetBinError(i+1));
        h_yrms_clY_24->SetBinError(i+1, h_y_clY_24->GetBinError(i+1)/TMath::Sqrt(h_y_clY_24->GetBinEntries(i+1)));

        h_xpullrms_clX_24->SetBinContent(i+1, h_pullx_clX_24->GetBinError(i+1));
        h_xpullrms_clX_24->SetBinError(i+1, h_pullx_clX_24->GetBinError(i+1)/TMath::Sqrt(h_pullx_clX_24->GetBinEntries(i+1)));

        h_ypullrms_clY_24->SetBinContent(i+1, h_pully_clY_24->GetBinError(i+1));
        h_ypullrms_clY_24->SetBinError(i+1, h_pully_clY_24->GetBinError(i+1)/TMath::Sqrt(h_pully_clY_24->GetBinEntries(i+1)));
    }
    
    // Now we do a lot of plotting.

    TCanvas *c1 = new TCanvas("c1", "", 0, 0, 800, 800);
    h_x->Draw("pe");
    h_x->Fit("gaus");
    h_x->GetFunction("gaus")->SetLineColor(kBlue);
    c1->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_residx.png");

    TCanvas *c2 = new TCanvas("c2", "", 0,0 , 800, 800);
    h_y->Draw("pe");
    h_y->Fit("gaus");
    h_y->GetFunction("gaus")->SetLineColor(kBlue);
    c2->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_residy.png");

    TCanvas *c1e = new TCanvas("c1e", "", 0, 0, 800, 800);
    h_xe->Draw("hist");
    // h_xe->Fit("gaus");
    // h_xe->GetFunction("gaus")->SetLineColor(kBlue);
    c1e->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_errorx.png");

    TCanvas *c2e = new TCanvas("c2e", "", 0, 0, 800, 800);
    h_ye->Draw("hist");
    // h_ye->Fit("gaus");
    // h_ye->GetFunction("gaus")->SetLineColor(kBlue);
    c2e->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_errory.png");

    TCanvas *c1_ClCheckX = new TCanvas("c1_ClCheckX", "", 0, 0, 800, 800);
    h_ClCheck_x->Draw("hist");
    // h_x->Fit("gaus");
    // h_x->GetFunction("gaus")->SetLineColor(kBlue);
    c1_ClCheckX->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_ClCheckX.png");

    TCanvas *c2_ClCheckY = new TCanvas("c2_ClCheckY", "", 0, 0, 800, 800);
    h_ClCheck_y->Draw("hist");
    // h_y->Fit("gaus");
    // h_y->GetFunction("gaus")->SetLineColor(kBlue);
    c2_ClCheckY->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_ClCheckY.png");

    TCanvas *c3 = new TCanvas("c3", "", 0, 0, 800, 800);
    h_pullx->Draw("pe");
    h_pullx->Fit("gaus");
    h_pullx->GetFunction("gaus")->SetLineColor(kBlue);
    c3->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_pullx.png");

    TCanvas *c4 = new TCanvas("c4", "", 0,0 , 800, 800);
    h_pully->Draw("pe");
    h_pully->Fit("gaus");
    h_pully->GetFunction("gaus")->SetLineColor(kBlue);
    c4->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_pully.png");

    TCanvas *c1_24 = new TCanvas("c1_24", "", 0, 0, 800, 800);
    h_x_24->Draw("pe");
    h_x_24->Fit("gaus");
    h_x_24->GetFunction("gaus")->SetLineColor(kBlue);
    c1_24->SaveAs("2024_layer1_fitChargePerSizeSkipPol_profile_residx.png");

    TCanvas *c2_24 = new TCanvas("c2_24", "", 0, 0, 800, 800);
    h_y_24->Draw("pe");
    h_y_24->Fit("gaus");
    h_y_24->GetFunction("gaus")->SetLineColor(kBlue);
    c2_24->SaveAs("2024_layer1_fitChargePerSizeSkipPol_profile_residy.png");

    TCanvas *c1e_24 = new TCanvas("c1e_24", "", 0, 0, 800, 800);
    h_xe_24->Draw("hist");
    // h_xe_24->Fit("gaus");
    // h_xe_24->GetFunction("gaus")->SetLineColor(kBlue);
    c1e_24->SaveAs("2024_layer1_fitChargePerSizeSkipPol_profile_errorx.png");

    TCanvas *c2e_24 = new TCanvas("c2e_24", "", 0, 0, 800, 800);
    h_ye_24->Draw("hist");
    // h_ye_24->Fit("gaus");
    // h_ye_24->GetFunction("gaus")->SetLineColor(kBlue);
    c2e_24->SaveAs("2024_layer1_fitChargePerSizeSkipPol_profile_errory.png");

    TCanvas *c1_ClCheckX_24 = new TCanvas("c1_ClCheckX_24", "", 0, 0, 800, 800);
    h_ClCheck_x_24->Draw("hist");
    // h_x->Fit("gaus");
    // h_x->GetFunction("gaus")->SetLineColor(kBlue);
    c1_ClCheckX_24->SaveAs("2024_layer1_fitChargePerSizeSkipPol_profile_ClCheckX.png");

    TCanvas *c2_ClCheckY_24 = new TCanvas("c2_ClCheckY_24", "", 0, 0, 800, 800);
    h_ClCheck_y_24->Draw("hist");
    // h_y->Fit("gaus");
    // h_y->GetFunction("gaus")->SetLineColor(kBlue);
    c2_ClCheckY_24->SaveAs("2024_layer1_fitChargePerSizeSkipPol_profile_ClCheckY.png");

    TCanvas *c3_24 = new TCanvas("c3_24", "", 0, 0, 800, 800);
    h_pullx_24->Draw("pe");
    h_pullx_24->Fit("gaus");
    h_pullx_24->GetFunction("gaus")->SetLineColor(kBlue);
    c3_24->SaveAs("2024_layer1_fitChargePerSizeSkipPol_profile_pullx.png");

    TCanvas *c4_24 = new TCanvas("c4_24", "", 0, 0, 800, 800);
    h_pully_24->Draw("pe");
    h_pully_24->Fit("gaus");
    h_pully_24->GetFunction("gaus")->SetLineColor(kBlue);
    c4_24->SaveAs("2024_layer1_fitChargePerSizeSkipPol_profile_pully.png");

    // TCanvas *c1_clx = new TCanvas("c1_clx", "", 0, 0, 800, 800);
    // h_x_clX->Draw();
    // h_x_clX->Fit("pol1");
    // h_x_clX->GetFunction("pol1")->SetLineColor(kBlue);
    // // h_x_clX->Fit("pol2","+");
    // // h_x_clX->GetFunction("pol2")->SetLineColor(kRed);
    // c1_clx->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_residx_clx.png");

    // TCanvas *c2_cly = new TCanvas("c2_cly", "", 0,0 , 800, 800);
    // h_y_clY->Draw("");
    // h_y_clY->Fit("pol1");
    // h_y_clY->GetFunction("pol1")->SetLineColor(kBlue);
    // // h_y_clY->Fit("pol2", "+");
    // // h_y_clY->GetFunction("pol2")->SetLineColor(kRed);
    // c2_cly->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_residy_cly.png");

    // TCanvas *c1_clxe = new TCanvas("c1_clxe", "", 0, 0, 800, 800);
    // h_xe_clX->SetMarkerStyle(20);
    // h_xe_clX->SetMinimum(0);
    // h_xe_clX->SetMaximum(100);
    // h_xe_clX->SetTitle("X Residual Error v. Cluster Size layer1");
    // h_xe_clX_24->SetMarkerStyle(20);
    // h_xe_clX_24->SetMarkerColor(kRed);
    // h_xe_clX->SetStats(0);
    // h_xe_clX_24->SetStats(0);
    // auto leg3 = new TLegend(0.7, 0.8, .9, .9);
    // leg3->AddEntry(h_xe_clX, "2021", "p");
    // leg3->AddEntry(h_xe_clX_24, "2024", "p");
    // h_xe_clX->Draw();
    // h_xe_clX_24->Draw("same");
    // leg3->Draw();
    // // h_xe_clX->Fit("pol1");
    // // h_xe_clX->GetFunction("pol1")->SetLineColor(kBlue);
    // // h_xe_clX->Fit("pol2", "+");
    // // h_xe_clX->GetFunction("pol2")->SetLineColor(kRed);
    // c1_clxe->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_residxe_clx.png");

    // TCanvas *c2_clye = new TCanvas("c2_clye", "", 0,0 , 800, 800);
    // h_ye_clY->SetMarkerStyle(20);
    // h_ye_clY->SetMinimum(0);
    // h_ye_clY->SetMaximum(100);
    // h_ye_clY->SetTitle("Y Residual Error v. Cluster Size layer1");
    // h_ye_clY_24->SetMarkerStyle(20);
    // h_ye_clY_24->SetMarkerColor(kRed);
    // h_ye_clY->SetStats(0);
    // h_ye_clY_24->SetStats(0);
    // auto leg2 = new TLegend(0.7, 0.8, .9, .9);
    // leg2->AddEntry(h_ye_clY, "2021", "p");
    // leg2->AddEntry(h_ye_clY_24, "2024", "p");
    // h_ye_clY->Draw();
    // h_ye_clY_24->Draw("same");
    // leg2->Draw();
    // // h_ye_clY->Fit("pol1");
    // // h_ye_clY->GetFunction("pol1")->SetLineColor(kBlue);
    // // h_ye_clY->Fit("pol2", "+");
    // // h_ye_clY->GetFunction("pol2")->SetLineColor(kRed);
    // c2_clye->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_residye_cly.png");

    // TCanvas *c3_clx = new TCanvas("c3_clx", "", 0, 0, 800, 800);
    // h_pullx_clX->Draw();
    // h_pullx_clX->Fit("pol1");
    // h_pullx_clX->GetFunction("pol1")->SetLineColor(kBlue);
    // // h_pullx_clX->Fit("pol2", "+");
    // // h_pullx_clX->GetFunction("pol2")->SetLineColor(kRed);
    // c3_clx->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_pullx_clx.png");

    // TCanvas *c4_cly = new TCanvas("c4_cly", "", 0,0 , 800, 800);
    // h_pully_clY->Draw();
    // h_pully_clY->Fit("pol1");
    // h_pully_clY->GetFunction("pol1")->SetLineColor(kBlue);
    // // h_pully_clY->Fit("pol2", "+");
    // // h_pully_clY->GetFunction("pol2")->SetLineColor(kRed);
    // c4_cly->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_pully_cly.png");

    // TCanvas *c1_clxe_rms = new TCanvas("c1_clxe_rms", "", 0, 0, 800, 800);
    // h_xpullrms_clX_qbin->SetMaximum(100);
    // h_xpullrms_clX_qbin->Draw("colz");
    // // h_x->Fit("gaus");
    // // h_x->GetFunction("gaus")->SetLineColor(kBlue);
    // c1_clxe_rms->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_pullxrms_clx_qbin.png");

    // TCanvas *c2_clye_rms = new TCanvas("c2_clye_rms", "", 0,0 , 800, 800);
    // h_ypullrms_clY_qbin->SetMaximum(20);
    // h_ypullrms_clY_qbin->Draw("colz");
    // // h_y->Fit("gaus");
    // // h_y->GetFunction("gaus")->SetLineColor(kBlue);
    // c2_clye_rms->SaveAs("2021_layer1_fitChargePerSizeSkipPol_profile_pullyrms_cly_qbin.png");

    // TCanvas *c1_clxrms = new TCanvas("c1_clxrms", "", 0, 0, 800, 800);
    // h_xrms_clX->SetMarkerStyle(20);
    // h_xrms_clX->SetMinimum(0);
    // h_xrms_clX->SetMaximum(50);
    // h_xrms_clX->SetTitle("X Residual RMS v. Cluster Size layer1");
    // h_xrms_clX_24->SetMarkerStyle(20);
    // h_xrms_clX_24->SetMarkerColor(kRed);
    // h_xrms_clX->Draw("pe");
    // h_xrms_clX_24->Draw("pe same");
    // // h_xrms_clX->Fit("pol2");
    // // h_xrms_clX->GetFunction("pol2")->SetLineColor(kBlue);
    // c1_clxrms->SaveAs("layer1_fitChargePerSizeSkipPol_profile_resxrms_clx.png");

    // TCanvas *c2_clyrms = new TCanvas("c2_clyrms", "", 0,0 , 800, 800);
    // h_yrms_clY->SetMarkerStyle(20);
    // h_yrms_clY->SetMinimum(0);
    // h_yrms_clY->SetMaximum(100);
    // h_yrms_clY->SetTitle("Y Residual RMS v. Cluster Size layer1");
    // h_yrms_clY_24->SetMarkerStyle(20);
    // h_yrms_clY_24->SetMarkerColor(kRed);
    // h_yrms_clY->Draw("pe");
    // h_yrms_clY_24->Draw("pe same");
    // // h_yrms_clY->Fit("pol2");
    // // h_yrms_clY->GetFunction("pol2")->SetLineColor(kBlue);
    // c2_clyrms->SaveAs("layer1_fitChargePerSizeSkipPol_profile_resyrms_cly.png");

    TCanvas *c1_clxpullrms = new TCanvas("c1_clxpullrms", "", 0, 0, 800, 800);
    h_xpullrms_clX->SetMarkerStyle(20);
    // h_xpullrms_clX->SetStats(0);
    h_xpullrms_clX->SetMinimum(0);
    h_xpullrms_clX->SetMaximum(5);
    h_xpullrms_clX->SetTitle("X Pull RMS v. Cluster Size layer1 2021");
    // h_xpullrms_clX_24->SetStats(0);
    // h_xpullrms_clX_24->SetMarkerStyle(20);
    // h_xpullrms_clX_24->SetMarkerColor(kRed);
    // auto leg = new TLegend(0.7, 0.8, .9, .9);
    // leg->AddEntry(h_xpullrms_clX, "2021", "p");
    // leg->AddEntry(h_xpullrms_clX_24, "2024", "p");
    h_xpullrms_clX->Draw("pe");
    // h_xpullrms_clX_24->Draw("pe same");
    // leg->Draw();
    // TF1 *f1 = new TF1("f1","pol2",0,10);
    // TF1 *f2 = new TF1("f2","pol2",0,10);
    h_xpullrms_clX->Fit("pol3");
    // h_xrms_clX->GetFunction("f1")->SetLineColor(kBlue);
    // h_xpullrms_clX_24->Fit("pol2","sames");
    // h_xrms_clX_24->GetFunction("f2")->SetLineColor(kRed);
    c1_clxpullrms->SaveAs("layer1_fitChargePerSizeSkipPol_profile_pullxrms_clx.png");

    TCanvas *c1_clxpullrms_24 = new TCanvas("c1_clxpullrms_24", "", 0, 0, 800, 800);
    // h_xpullrms_clX->SetStats(0);
    h_xpullrms_clX_24->SetMinimum(0);
    h_xpullrms_clX_24->SetMaximum(5);
    h_xpullrms_clX_24->SetTitle("X Pull RMS v. Cluster Size layer1 2024");
    // h_xpullrms_clX_24->SetStats(0);
    h_xpullrms_clX_24->SetMarkerStyle(20);
    h_xpullrms_clX_24->SetMarkerColor(kRed);
    // auto leg = new TLegend(0.7, 0.8, .9, .9);
    // leg->AddEntry(h_xpullrms_clX, "2021", "p");
    // leg->AddEntry(h_xpullrms_clX_24, "2024", "p");
    h_xpullrms_clX_24->Draw("pe same");
    // leg->Draw();
    // TF1 *f1 = new TF1("f1","pol2",0,10);
    // TF1 *f2 = new TF1("f2","pol2",0,10);
    // h_xrms_clX->GetFunction("f1")->SetLineColor(kBlue);
    h_xpullrms_clX_24->Fit("pol2", "sames");
    // h_xrms_clX_24->GetFunction("f2")->SetLineColor(kRed);
    c1_clxpullrms_24->SaveAs("layer1_fitChargePerSizeSkipPol_profile_pullxrms_clx_24.png");

    TCanvas *c2_clypullrms = new TCanvas("c2_clypullrms", "", 0, 0, 800, 800);
    h_ypullrms_clY->SetMarkerStyle(20);
    // h_ypullrms_clY->SetStats(0);
    h_ypullrms_clY->SetMinimum(0);
    h_ypullrms_clY->SetMaximum(5);
    h_ypullrms_clY->SetTitle("Y Pull RMS v. Cluster Size layer1 2021");
    // h_ypullrms_clY_24->SetMarkerStyle(20);
    // h_ypullrms_clY_24->SetStats(0);
    // h_ypullrms_clY_24->SetMarkerColor(kRed);
    // auto leg1 = new TLegend(0.7, 0.8, .9, .9);
    // leg1->AddEntry(h_ypullrms_clY, "2021", "p");
    // leg1->AddEntry(h_ypullrms_clY_24, "2024", "p");
    h_ypullrms_clY->Draw("pe ");
    // h_ypullrms_clY_24->Draw("pe same");
    // leg1->Draw();
    h_ypullrms_clY->Fit("pol2");
    // h_yrms_clY->GetFunction("pol2")->SetLineColor(kBlue);
    // h_ypullrms_clY_24->Fit("pol2","sames");
    // h_yrms_clY_24->GetFunction("pol2")->SetLineColor(kRed);
    c2_clypullrms->SaveAs("layer1_fitChargePerSizeSkipPol_profile_pullyrms_cly.png");

    TCanvas *c2_clypullrms_24 = new TCanvas("c2_clypullrms_24", "", 0, 0, 800, 800);
    // h_ypullrms_clY->SetMarkerStyle(20);
    // h_ypullrms_clY->SetStats(0);
    h_ypullrms_clY_24->SetMinimum(0);
    h_ypullrms_clY_24->SetMaximum(5);
    h_ypullrms_clY_24->SetTitle("Y Pull RMS v. Cluster Size layer1 2024");
    h_ypullrms_clY_24->SetMarkerStyle(20);
    // h_ypullrms_clY_24->SetStats(0);
    h_ypullrms_clY_24->SetMarkerColor(kRed);
    // auto leg1 = new TLegend(0.7, 0.8, .9, .9);
    // leg1->AddEntry(h_ypullrms_clY, "2021", "p");
    // leg1->AddEntry(h_ypullrms_clY_24, "2024", "p");
    // h_ypullrms_clY->Draw("pe ");
    h_ypullrms_clY_24->Draw("pe same");
    // leg1->Draw();
    // h_ypullrms_clY->Fit("pol3");
    // h_yrms_clY->GetFunction("pol2")->SetLineColor(kBlue);
    h_ypullrms_clY_24->Fit("pol2", "sames");
    // h_yrms_clY_24->GetFunction("pol2")->SetLineColor(kRed);
    c2_clypullrms_24->SaveAs("layer1_fitChargePerSizeSkipPol_profile_pullyrms_cly_24.png");

    h_x->Reset(); h_y->Reset(); h_pullx->Reset(); h_pully->Reset(); h_x_clX->Reset(); h_y_clY->Reset(); h_xe_clX->Reset(); h_ye_clY->Reset(); h_pullx_clX->Reset(); h_pully_clY->Reset();
    h_xrms_clX->Reset(); h_yrms_clY->Reset(); h_xpullrms_clX_qbin->Reset(); h_ypullrms_clY_qbin->Reset();

    h_x_24->Reset();
    h_y_24->Reset();
    h_pullx_24->Reset();
    h_pully_24->Reset();
    h_x_clX_24->Reset();
    h_y_clY_24->Reset();
    h_xe_clX_24->Reset();
    h_ye_clY_24->Reset();
    h_pullx_clX_24->Reset();
    h_pully_clY_24->Reset();
    h_xrms_clX_24->Reset();
    h_yrms_clY_24->Reset();
    h_xpullrms_clX_qbin_24->Reset();
    h_ypullrms_clY_qbin_24->Reset();
    //Edge clusters

    // use_edge = true;
    // fill_pulls(f_name, h_x, h_y, h_pullx, h_pully, use_edge, use_badpix);

    // TCanvas *c1e = new TCanvas("c1", "", 0, 0, 800, 800);
    // h_x->Draw("pe");
    // h_x->Fit("gaus");
    // h_x->GetFunction("gaus")->SetLineColor(kBlue);
    // c1e->SaveAs("edge_residx.png");

    // TCanvas *c2e = new TCanvas("c2", "", 0,0 , 800, 800);
    // h_y->Draw("pe");
    // h_y->Fit("gaus");
    // h_y->GetFunction("gaus")->SetLineColor(kBlue);
    // c2e->SaveAs("edge_residy.png");

    // TCanvas *c3e = new TCanvas("c3", "", 0, 0, 800, 800);
    // h_pullx->Draw("pe");
    // h_pullx->Fit("gaus");
    // h_pullx->GetFunction("gaus")->SetLineColor(kBlue);
    // c3e->SaveAs("edge_pullx.png");

    // TCanvas *c4e = new TCanvas("c4", "", 0,0 , 800, 800);
    // h_pully->Draw("pe");
    // h_pully->Fit("gaus");
    // h_pully->GetFunction("gaus")->SetLineColor(kBlue);
    // c4e->SaveAs("edge_pully.png");

    // h_x->Reset(); h_y->Reset(); h_pullx->Reset(); h_pully->Reset();
    // //BadPix clusters
    // use_edge = false;
    // use_badpix = true;
    // fill_pulls(f_name, h_x, h_y, h_pullx, h_pully, use_edge, use_badpix);

    // TCanvas *c1b = new TCanvas("c1", "", 0, 0, 800, 800);
    // h_x->Draw("pe");
    // h_x->Fit("gaus");
    // h_x->GetFunction("gaus")->SetLineColor(kBlue);
    // c1b->SaveAs("badpix_residx.png");

    // TCanvas *c2b = new TCanvas("c2", "", 0,0 , 800, 800);
    // h_y->Draw("pe");
    // h_y->Fit("gaus");
    // h_y->GetFunction("gaus")->SetLineColor(kBlue);
    // c2b->SaveAs("badpix_residy.png");

    // TCanvas *c3b = new TCanvas("c3", "", 0, 0, 800, 800);
    // h_pullx->Draw("pe");
    // h_pullx->Fit("gaus");
    // h_pullx->GetFunction("gaus")->SetLineColor(kBlue);
    // c3b->SaveAs("badpix_pullx.png");

    // TCanvas *c4b = new TCanvas("c4", "", 0,0 , 800, 800);
    // h_pully->Draw("pe");
    // h_pully->Fit("gaus");
    // h_pully->GetFunction("gaus")->SetLineColor(kBlue);
    // c4b->SaveAs("badpix_pully.png");

    return;
}