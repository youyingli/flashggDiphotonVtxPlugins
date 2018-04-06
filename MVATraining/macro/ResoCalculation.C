
double getReso( TH1* hist, double critiria = 0.683) {

    int peak = 0;
    int peakbin = hist->FindBin(peak);

    int minus_bin = -1;
    int plus_bin = -1;
    for (int i = 1; i <= hist->GetNbinsX() / 2; i++) {
        double area = hist->Integral(peakbin - i, peakbin + i);
        double ratio_area = area / hist->GetEntries();
        if (ratio_area > critiria) {
            minus_bin = peakbin - i;
            plus_bin = peakbin + i;
            break;
        }
    }
    double minus_val = fabs(hist->GetXaxis()->GetBinCenter(minus_bin));
    double plus_val = hist->GetXaxis()->GetBinCenter(plus_bin);
    return 0.5 * (minus_val + plus_val);

}

void ResoCalculation() {

    int ntracks;   
    int nphotons;   
    int photon_index;   
    double dipho_mass;   
    double pho_pt;   
    double ztrue    ;
    double zvtx     ;
    double z1Pix    ;
    double z1Tib    ;
    double z1Tob    ;
    double z1PixFwd ;
    double z1Tid    ;
    double z1Tec    ;
    double z2Pix    ;
    double z2Tib    ;
    double z2Tob    ;
    double z2PixFwd ;
    double z2Tid    ;
    double z2Tec    ;
    
    double factor = 5;
    double nbin = 10000;
    double sigma[24] = {0};
    sigma[0]   = 0.0125255;
    sigma[1]   = 0.716301;
    sigma[2]   = 3.17615;
    sigma[3]   = 0.0581667;
    sigma[4]   = 0.38521;
    sigma[5]   = 1.67937;
    sigma[6]   = 0.0298574;
    sigma[7]   = 0.414393;
    sigma[8]   = 1.06805;
    sigma[9]   = 0.180419;
    sigma[10]   = 0.494722;
    sigma[11]   = 1.21941;
    sigma[12]   = 0.0178107;
    sigma[13]   = 1.3188;
    sigma[14]   = 2.23662;
    sigma[15]   = 0.152157;
    sigma[16]   = 0.702755;
    sigma[17]   = 2.46599;
    sigma[18]   = 0.0935307;
    sigma[19]   = 0.756568;
    sigma[20]   = 0.62143;
    sigma[21]   = 0.577081;
    sigma[22]   = 0.892751;
    sigma[23]   = 1.56638;

    TFile* file = TFile::Open("resolution.root");
    TTree* tree = (TTree*) file->Get("commissioning/resolutionTree");
    tree->SetBranchAddress("ntracks"     , &ntracks     );
    tree->SetBranchAddress("nphotons"    , &nphotons    );
    tree->SetBranchAddress("photon_index", &photon_index);
    tree->SetBranchAddress("dipho_mass"  , &dipho_mass  );
    tree->SetBranchAddress("pho_pt"      , &pho_pt      );
    tree->SetBranchAddress("ztrue"       , &ztrue       );
    tree->SetBranchAddress("zvtx"        , &zvtx        );
    tree->SetBranchAddress("z1Pix"       , &z1Pix       );
    tree->SetBranchAddress("z1Tib"       , &z1Tib       );
    tree->SetBranchAddress("z1Tob"       , &z1Tob       );
    tree->SetBranchAddress("z1PixFwd"    , &z1PixFwd    );
    tree->SetBranchAddress("z1Tid"       , &z1Tid       );
    tree->SetBranchAddress("z1Tec"       , &z1Tec       );
    tree->SetBranchAddress("z2Pix"       , &z2Pix       );
    tree->SetBranchAddress("z2Tib"       , &z2Tib       );
    tree->SetBranchAddress("z2Tob"       , &z2Tob       );
    tree->SetBranchAddress("z2PixFwd"    , &z2PixFwd    );
    tree->SetBranchAddress("z2Tid"       , &z2Tid       );
    tree->SetBranchAddress("z2Tec"       , &z2Tec       );

    TFile* outfile = TFile::Open("resout.root", "recreate");

    TH1F* h1_sigma1Pix     = new TH1F("h1_sigma1Pix",    "", nbin, -1 * factor * sigma[0], factor * sigma[0]);
    TH1F* h1_sigma1Tib     = new TH1F("h1_sigma1Tib",    "", nbin, -1 * factor * sigma[1], factor * sigma[1]);
    TH1F* h1_sigma1Tob     = new TH1F("h1_sigma1Tob",    "", nbin, -1 * factor * sigma[2], factor * sigma[2]);
    TH1F* h1_sigma1PixFwd  = new TH1F("h1_sigma1PixFwd", "", nbin, -1 * factor * sigma[3], factor * sigma[3]);
    TH1F* h1_sigma1Tid     = new TH1F("h1_sigma1Tid",    "", nbin, -1 * factor * sigma[4], factor * sigma[4]);
    TH1F* h1_sigma1Tec     = new TH1F("h1_sigma1Tec",    "", nbin, -1 * factor * sigma[5], factor * sigma[5]);
    TH1F* h1_sigma2Pix     = new TH1F("h1_sigma2Pix",    "", nbin, -1 * factor * sigma[6], factor * sigma[6]);
    TH1F* h1_sigma2Tib     = new TH1F("h1_sigma2Tib",    "", nbin, -1 * factor * sigma[7], factor * sigma[7]);
    TH1F* h1_sigma2Tob     = new TH1F("h1_sigma2Tob",    "", nbin, -1 * factor * sigma[8], factor * sigma[8]);
    TH1F* h1_sigma2PixFwd  = new TH1F("h1_sigma2PixFwd", "", nbin, -1 * factor * sigma[9], factor * sigma[9]);
    TH1F* h1_sigma2Tid     = new TH1F("h1_sigma2Tid",    "", nbin, -1 * factor * sigma[10], factor * sigma[10]);
    TH1F* h1_sigma2Tec     = new TH1F("h1_sigma2Tec",    "", nbin, -1 * factor * sigma[11], factor * sigma[11]);

    TH1F* h1_singlelegsigma1Pix     = new TH1F("h1_singlelegsigma1Pix",    "", nbin, -1 * factor * sigma[12], factor * sigma[12]);
    TH1F* h1_singlelegsigma1Tib     = new TH1F("h1_singlelegsigma1Tib",    "", nbin, -1 * factor * sigma[13], factor * sigma[13]);
    TH1F* h1_singlelegsigma1Tob     = new TH1F("h1_singlelegsigma1Tob",    "", nbin, -1 * factor * sigma[14], factor * sigma[14]);
    TH1F* h1_singlelegsigma1PixFwd  = new TH1F("h1_singlelegsigma1PixFwd", "", nbin, -1 * factor * sigma[15], factor * sigma[15]);
    TH1F* h1_singlelegsigma1Tid     = new TH1F("h1_singlelegsigma1Tid",    "", nbin, -1 * factor * sigma[16], factor * sigma[16]);
    TH1F* h1_singlelegsigma1Tec     = new TH1F("h1_singlelegsigma1Tec",    "", nbin, -1 * factor * sigma[17], factor * sigma[17]);
    TH1F* h1_singlelegsigma2Pix     = new TH1F("h1_singlelegsigma2Pix",    "", nbin, -1 * factor * sigma[18], factor * sigma[18]);
    TH1F* h1_singlelegsigma2Tib     = new TH1F("h1_singlelegsigma2Tib",    "", nbin, -1 * factor * sigma[19], factor * sigma[19]);
    TH1F* h1_singlelegsigma2Tob     = new TH1F("h1_singlelegsigma2Tob",    "", nbin, -1 * factor * sigma[20], factor * sigma[20]);
    TH1F* h1_singlelegsigma2PixFwd  = new TH1F("h1_singlelegsigma2PixFwd", "", nbin, -1 * factor * sigma[21], factor * sigma[21]);
    TH1F* h1_singlelegsigma2Tid     = new TH1F("h1_singlelegsigma2Tid",    "", nbin, -1 * factor * sigma[22], factor * sigma[22]);
    TH1F* h1_singlelegsigma2Tec     = new TH1F("h1_singlelegsigma2Tec",    "", nbin, -1 * factor * sigma[23], factor * sigma[23]);

    double minval = -500.;
    double maxval = 500.;
    for (int ientry = 0; ientry < tree->GetEntries(); ientry++) {
        tree->GetEntry(ientry);
        if (nphotons != 2) continue;

        if (ntracks == 2) {

            if ( z1Pix    > minval && z1Pix    < maxval )   h1_sigma1Pix    ->Fill(zvtx - z1Pix    ); 
            if ( z1Tib    > minval && z1Tib    < maxval )   h1_sigma1Tib    ->Fill(zvtx - z1Tib    ); 
            if ( z1Tob    > minval && z1Tob    < maxval )   h1_sigma1Tob    ->Fill(zvtx - z1Tob    ); 
            if ( z1PixFwd > minval && z1PixFwd < maxval )   h1_sigma1PixFwd ->Fill(zvtx - z1PixFwd ); 
            if ( z1Tid    > minval && z1Tid    < maxval )   h1_sigma1Tid    ->Fill(zvtx - z1Tid    ); 
            if ( z1Tec    > minval && z1Tec    < maxval )   h1_sigma1Tec    ->Fill(zvtx - z1Tec    ); 
            if ( z2Pix    > minval && z2Pix    < maxval )   h1_sigma2Pix    ->Fill(zvtx - z2Pix    ); 
            if ( z2Tib    > minval && z2Tib    < maxval )   h1_sigma2Tib    ->Fill(zvtx - z2Tib    ); 
            if ( z2Tob    > minval && z2Tob    < maxval )   h1_sigma2Tob    ->Fill(zvtx - z2Tob    ); 
            if ( z2PixFwd > minval && z2PixFwd < maxval )   h1_sigma2PixFwd ->Fill(zvtx - z2PixFwd ); 
            if ( z2Tid    > minval && z2Tid    < maxval )   h1_sigma2Tid    ->Fill(zvtx - z2Tid    ); 
            if ( z2Tec    > minval && z2Tec    < maxval )   h1_sigma2Tec    ->Fill(zvtx - z2Tec    ); 

        } else if (ntracks == 1) {

            if ( z1Pix    > minval && z1Pix    < maxval )  h1_singlelegsigma1Pix    ->Fill(zvtx - z1Pix    ); 
            if ( z1Tib    > minval && z1Tib    < maxval )  h1_singlelegsigma1Tib    ->Fill(zvtx - z1Tib    ); 
            if ( z1Tob    > minval && z1Tob    < maxval )  h1_singlelegsigma1Tob    ->Fill(zvtx - z1Tob    ); 
            if ( z1PixFwd > minval && z1PixFwd < maxval )  h1_singlelegsigma1PixFwd ->Fill(zvtx - z1PixFwd ); 
            if ( z1Tid    > minval && z1Tid    < maxval )  h1_singlelegsigma1Tid    ->Fill(zvtx - z1Tid    ); 
            if ( z1Tec    > minval && z1Tec    < maxval )  h1_singlelegsigma1Tec    ->Fill(zvtx - z1Tec    ); 
            if ( z2Pix    > minval && z2Pix    < maxval )  h1_singlelegsigma2Pix    ->Fill(zvtx - z2Pix    ); 
            if ( z2Tib    > minval && z2Tib    < maxval )  h1_singlelegsigma2Tib    ->Fill(zvtx - z2Tib    ); 
            if ( z2Tob    > minval && z2Tob    < maxval )  h1_singlelegsigma2Tob    ->Fill(zvtx - z2Tob    ); 
            if ( z2PixFwd > minval && z2PixFwd < maxval )  h1_singlelegsigma2PixFwd ->Fill(zvtx - z2PixFwd ); 
            if ( z2Tid    > minval && z2Tid    < maxval )  h1_singlelegsigma2Tid    ->Fill(zvtx - z2Tid    ); 
            if ( z2Tec    > minval && z2Tec    < maxval )  h1_singlelegsigma2Tec    ->Fill(zvtx - z2Tec    ); 

        }

    }

    cout << getReso( h1_sigma1Pix   ) << endl; 
    cout << getReso( h1_sigma1Tib   ) << endl; 
    cout << getReso( h1_sigma1Tob   ) << endl; 
    cout << getReso( h1_sigma1PixFwd) << endl; 
    cout << getReso( h1_sigma1Tid   ) << endl; 
    cout << getReso( h1_sigma1Tec   ) << endl; 
    cout << getReso( h1_sigma2Pix   ) << endl; 
    cout << getReso( h1_sigma2Tib   ) << endl; 
    cout << getReso( h1_sigma2Tob   ) << endl; 
    cout << getReso( h1_sigma2PixFwd) << endl; 
    cout << getReso( h1_sigma2Tid   ) << endl; 
    cout << getReso( h1_sigma2Tec   ) << endl; 

    cout << getReso( h1_singlelegsigma1Pix   ) << endl; 
    cout << getReso( h1_singlelegsigma1Tib   ) << endl; 
    cout << getReso( h1_singlelegsigma1Tob   ) << endl; 
    cout << getReso( h1_singlelegsigma1PixFwd) << endl; 
    cout << getReso( h1_singlelegsigma1Tid   ) << endl; 
    cout << getReso( h1_singlelegsigma1Tec   ) << endl; 
    cout << getReso( h1_singlelegsigma2Pix   ) << endl; 
    cout << getReso( h1_singlelegsigma2Tib   ) << endl; 
    cout << getReso( h1_singlelegsigma2Tob   ) << endl; 
    cout << getReso( h1_singlelegsigma2PixFwd) << endl; 
    cout << getReso( h1_singlelegsigma2Tid   ) << endl; 
    cout << getReso( h1_singlelegsigma2Tec   ) << endl; 

    outfile->Write();
    outfile->Close();
    file->Close();
    
}
