//gamma n -> K^+ Sigma^-; Sigma^- p -> Lambda n
{
  gStyle->SetOptStat(0);
  TFile fileOutput1("Sigma_Proton_Output_test.root","recreate");
  // Creating histograms
  ////////////////////////////// Angular dependence histograms //////////////////////////////
  TH2D* h_Kaon_Plus=new TH2D("h_Kaon_Plus","Angular dependence of K^{+}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,6,200,0,180);
  TH2D* h_Sigma_Minus=new TH2D("h_Sigma_Minus","Angular dependence of #Sigma^{-}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,6,200,0,180);
  TH2D* h_Lambda=new TH2D("h_Lambda","Angular dependence of #Lambda^{0}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,6,200,0,180);
  TH2D* h_Neutron=new TH2D("h_Sigma","Angular dependence of n; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,6,200,0,180);
  TH2D* h_Proton=new TH2D("h_Proton","Angular dependence of Proton produced from #Lambda^{0} decay; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,6,200,0,180);  
  TH2D* h_Pi_Minus=new TH2D("h_Pi_Minus","Angular dependence of #pi^{-} produced from #Lambda^{0} decay; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,6,200,0,180);

  ////////////////////////////// Constrained Angular dependence histograms //////////////////////////////
  TH2D* h_cKaon_Plus=new TH2D("h_cKaon_Plus","Angular dependence of K^{+}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,6,200,0,180);
  TH2D* h_cSigma_Minus=new TH2D("h_cSigma_Minus","Angular dependence of #Sigma^{-}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,6,200,0,180);
  TH2D* h_cLambda=new TH2D("h_cLambda","Angular dependence of #Lambda^{0}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,6,200,0,180);
  TH2D* h_cNeutron=new TH2D("h_cSigma","Angular dependence of n; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,6,200,0,180);
  TH2D* h_cProton=new TH2D("h_cProton","Angular dependence of Proton produced from #Lambda^{0} decay; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,6,200,0,180);  
  TH2D* h_cPi_Minus=new TH2D("h_cPi_Minus","Angular dependence of #pi^{-} produced from #Lambda^{0} decay; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,6,200,0,180);

  ////////////////////////////// Other Histograms //////////////////////////////
  TH1F* h_Sigma_Momentum1 = new TH1F("h_Sigma_Momentum1","Momentum of #Sigma^{-}; Momentum [GeV/c]; Counts",12000,0,5);
  TH1F* h_Sigma_Momentum2 = new TH1F("h_Sigma_Momentum2","Momentum of #Sigma^{-}; Momentum [GeV/c]; Counts",12000,0,5);
  TH1F* h_Sigma_Scattered = new TH1F("h_Sigma_Scattered","Number of #Sigma^{-} Scattering Events per 100 days per MeV; Momentum [MeV/c]; Number of Events",12000,0,5);
  TH1F* h_Sigma_Momentum_Constrained = new TH1F("h_Sigma_Momentum_Constrained","Momentum of #Sigma^{-}; Momentum [GeV/c]; Counts",12000,0,5);
  TH1F* h_Sigma_Momentum = new TH1F("h_Sigma_Momentum","Momentum Distribution of #Sigma^{-}; Momentum [GeV/c]; Counts",12000,0,5); // Sigma_Minus->Rho()
  TH1F* h_Sigma_Acceptance = new TH1F("h_Sigma_Acceptance","Acceptance of #Sigma^{-}; Momentum [GeV/c]; Angle [#theta, degrees]",12000,0,5); // Acceptance as a function of Sigma_Momentum
  TH1F* h_Sigma_Detected = new TH1F("h_Sigma_Detected","Number of #Sigma^{-} detected over 100 days; Momentum [GeV/c]; Counts",12000,0,5); // Integral = Total number of sigma detected

  TH1F* h_Beam_Momentum1 = new TH1F("h_Beam_Momentum1","Momentum of beam w.r.t flux; Momentum [GeV/c]; Events per second",12000,0,12);
  TH1F* h_Beam_Momentum2 = new TH1F("h_Beam_Momentum2","Momentum of beam w.r.t flux; Momentum [GeV/c]; Events per second",12000,0,12);
  TH1F* h_Beam_Momentum3 = new TH1F("h_Beam_Momentum3","Momentum of beam w.r.t flux; Momentum [GeV/c]; Events per second",12000,0,12);
  TH1F* h_Beam_Flux = new TH1F("h_Beam_Flux","Photon Flux w.r.t Momentum; Momentum [GeV/c]; Photons per second",12000,0,12); // Beam Flux, no. photons/mev/s
  TH1F* h_Sigma_Cross_Section = new TH1F("h_Sigma_Cross_Section","Cross Section Distribution for Sigma Production; Momentum [GeV/c]; Cross Section [#mu b];",12000,0,12); // Sigma production cross section
  TH1F* h_Beam_Momentum6 = new TH1F("h_Beam_Momentum6","Momentum of beam w.r.t flux; Momentum [GeV/c]; Events per second",12000,0,12); // Sigma flux no. sigma/mev/s

  ////////////////////////////// Path Length Plots //////////////////////////////
  TF1* Decay_Plot = new TF1("Decay_Plot","TMath::Exp(-x/([0]*[1]*30*0.1479))", 0, 100);  // 0.1479 = Sigma^- mean lifetime in ns
  TH2D* h_l_xy = new TH2D("h_l_xy","#phi of #Sigma^{-} vs X-Y component of path length; #phi [radians]; length [cm]",200,-4,4,200,0,10);
  TH2D* h_l_z = new TH2D("h_l_z","#theta of #Sigma^{-} vs Z component of path length; #theta [degrees]; length [cm]",200,0,180,200,0,50);
  TH2D* h_Path_Length = new TH2D("h_Path_Length","#theta of #Sigma^{-} vs path length; #theta [degrees]; length [cm]",200,0,180,200,0,50);
  TH1F* h_Filled_Path_Length = new TH1F("h_Filled_Path_Length","Path length; Path Length [cm]; Counts ", 5000,0,40);
  TH2D* h_Path_Length_Decay = new TH2D("h_Path_Length_Decay","#theta of #Sigma^{-} decay vs path length; #theta [degrees]; length [cm]",200,0,180,200,0,50);
//------------------------------------------------------------------------------------------------------------------------------------------------------
  // Creating particles and terms
  // How many events to simulate and percentage completed
  Int_t nevents=1000000;
  Int_t Percentage=nevents/100;

  // Creating TLorentzVectors of particles
  TLorentzVector Beam, Target_Neutron, Target_Proton;       // Beam and target
  TLorentzVector *Kaon_Plus, *Sigma_Minus;                  // First vertex particles
  TLorentzVector *Lambda, *Neutron;                         // Second vertex particles
  TLorentzVector *Proton, *Pi_Minus;                        // Third Vertex Particles

  // Making Weights
  Double_t Phasespace_Weight_1, Phasespace_Weight_2, Phasespace_Weight_3;

  // Setting TLorentzVectors for beam and target in GeV (Px,Py,Pz,M)
  Target_Neutron.SetXYZM(0,0,0,0.93957);                    // M is 0.93957 (Neutron Target)
  Target_Proton.SetXYZM(0,0,0,0.93827);                     // M is 0.93827 (Proton Target)

  // Defining initial vertex and masses of particles
  TLorentzVector V1 = Beam + Target_Neutron;                // V1 = total energy, Beam is Gamma photons, Target is Neutron
  Double_t Masses_1[2] = {1.119745,0.49368};                // Sigma^-, K^+ (primary vertex)
  Double_t Masses_2[2] = {1.11568,0.93957};                 // Lambda n (Sigma^- interaction with proton)
  Double_t Masses_3[2] = {0.93827,0.13957};                 // Proton and Pi^- (Lambda Decay)


  // Creating decay vertices
  TGenPhaseSpace Vertex_1, Vertex_2, Vertex_3;

  TF1 *BeamEnergy = new TF1("BeamEnergy", "1.0/x", 0.9,5.0);

  ////////////////////////////////// Defining Variables For Path Length Calculations //////////////////////////////////
  // Target
  Double_t Target_Length = 40;                              // 40 cm
  Double_t Target_Radius = 3;                               // 3 cm
  // Path Length Components
  Double_t Path_Length;                                     // Path Length
  Double_t l_xy;                                            // Path Length x-y component
  Double_t l_z;                                             // Path Length z component
  // Angle Components
  Double_t Phi_v;
  Double_t Phi_Sigma;
  Double_t Theta_Sigma;
  Double_t Alpha;                                           
  // Vertex Components
  Double_t x_v;
  Double_t y_v;
  Double_t z_v;
  Double_t r_v;                                             // Radius from centre of target to vertex position
  Double_t z_exit;                                          // z vertex of exit position
  // Decay Components
  Double_t Sigma_Beta;
  Double_t Sigma_Gamma;
  Double_t Decay_Path;
  Double_t Path_Length_Decay;
  //--------------------------------------------------------------------------------------------------------------------------------------

  TRandom3 *myrndm = new TRandom3();
  
  for (Int_t i1=1;i1<12001;i1++){
    if (h_Beam_Flux->GetBinCenter(i1)>0.9 && h_Beam_Flux->GetBinCenter(i1) < 5)h_Beam_Flux->SetBinContent(i1,BeamEnergy->Eval(h_Beam_Flux->GetBinCenter(i1)));
    if (h_Sigma_Cross_Section->GetBinCenter(i1)>1.1 && h_Sigma_Cross_Section->GetBinCenter(i1) < 5)h_Sigma_Cross_Section->SetBinContent(i1,1.0);
  }

  Double_t Scale = 12000000.0/h_Beam_Flux->Integral();
  h_Beam_Flux->Scale(Scale);
  h_Beam_Momentum6->Multiply(h_Beam_Flux,h_Sigma_Cross_Section);
  h_Beam_Momentum6->Scale(4*0.16*3/1000000);

////////////////////////////////////////////////////////// Looping over simulated events  //////////////////////////////////////////////////////////
for (Long64_t i=0;i<nevents;i++){

  // Counter, shows percentage of simulated events completed
  if (i % Percentage == 0){
    fprintf (stderr, "%lld\r", i/Percentage);
    fflush (stderr);
  }

// Setting Decays --------------------------------------------------------------------------------------------------------------------------------------
  Beam.SetXYZM(0,0,BeamEnergy->GetRandom(),0);            // Pz is randomized, M is 0 (Photon Beam)
  V1 = Beam + Target_Neutron; 
  // .SetDecay(total energy, no. product particles, mass array)
  // Setting the first decay
  if (!Vertex_1.SetDecay(V1,2,Masses_1)) continue;
  // Generating event and defining the phasespace weight for this decay
  Phasespace_Weight_1=Vertex_1.Generate();
  // Assigning the decay particles from the array above
  Kaon_Plus=Vertex_1.GetDecay(0);
  Sigma_Minus=Vertex_1.GetDecay(1);

  // Getting total energy for second decay from Sigma^-
  TLorentzVector V2 = (TLorentzVector)*Sigma_Minus + Target_Proton;
  // Setting the second decay (Lambda n)
  Vertex_2.SetDecay(V2,2,Masses_2);
  // Defining the phasespace for this decay
  Phasespace_Weight_2 = Vertex_2.Generate();
  Lambda=Vertex_2.GetDecay(0);
  Neutron=Vertex_2.GetDecay(1);

  // Getting total energy for third decay from Lambda
  TLorentzVector V3 = (TLorentzVector)*Lambda;
  // Setting the third decay (Lambda decay to Proton + Pi^-)
  Vertex_3.SetDecay(V3,2,Masses_3);
  // Defining the phasespace for this decay
  Phasespace_Weight_3 = Vertex_3.Generate();
  Proton = Vertex_3.GetDecay(0);
  Pi_Minus = Vertex_3.GetDecay(1);
//------------------------------------------------------------------------------------------------------------------------------------------------------
  // Path Length Calculations
  r_v=5;
  while(r_v>3){
    x_v = myrndm->Gaus(0,1); // In cm
    y_v = myrndm->Gaus(0,1); // In cm
    r_v=sqrt(x_v*x_v+y_v*y_v);
  }
  Phi_v = TMath::ATan2(y_v,x_v);
  Phi_Sigma = Sigma_Minus->Phi();
  Theta_Sigma = Sigma_Minus->Theta();
  z_v = myrndm->Rndm()*Target_Length;
  Alpha = (TMath::ATan(1)*4) - Phi_v + Phi_Sigma;

  l_xy = (r_v*TMath::Cos(Alpha)) + sqrt((r_v*r_v*(TMath::Cos(Alpha))*(TMath::Cos(Alpha))) + Target_Radius*Target_Radius - r_v*r_v);
  l_z = l_xy/TMath::Tan(Theta_Sigma);
  z_exit = l_z + z_v;

  if (z_exit > Target_Length){
    l_z = Target_Length - z_v;
    l_xy = l_z*TMath::Tan(Theta_Sigma);
    Path_Length = l_z/TMath::Cos(Theta_Sigma);
  }
  else {
    Path_Length = l_xy/TMath::Sin(Theta_Sigma);
  }
  // Setting the decay constants for Sigma
  Sigma_Beta = Sigma_Minus->Beta();
  Sigma_Gamma = Sigma_Minus->Gamma();

  Decay_Plot->SetParameter(0, Sigma_Beta);
  Decay_Plot->SetParameter(1, Sigma_Gamma);
  Decay_Path = Decay_Plot->GetRandom();
  // Plotting path length decay lengths
  if (Decay_Path < Path_Length){
    Path_Length_Decay = Decay_Path;
  }
  else {
    Path_Length_Decay = Path_Length;
  }
  h_Filled_Path_Length->Fill(Path_Length_Decay);
//------------------------------------------------------------------------------------------------------------------------------------------------------
  h_Beam_Momentum1->Fill(Beam.Rho());
  h_Sigma_Momentum1->Fill(Sigma_Minus->Rho());
  h_Sigma_Momentum->Fill(Sigma_Minus->Rho());

  // Filling angular distribution histograms (x,y,weights applied)
  //Vertex 1
  h_Kaon_Plus->Fill(Kaon_Plus->Rho(),Kaon_Plus->Theta()*TMath::RadToDeg(),Phasespace_Weight_1);
  h_Sigma_Minus->Fill(Sigma_Minus->Rho(),Sigma_Minus->Theta()*TMath::RadToDeg(),Phasespace_Weight_1);
  //Vertex 2
  h_Lambda->Fill(Lambda->Rho(),Lambda->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2);
  h_Neutron->Fill(Neutron->Rho(),Neutron->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2);
  //Vertex 3
  h_Proton->Fill(Proton->Rho(),Proton->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3);
  h_Pi_Minus->Fill(Pi_Minus->Rho(),Pi_Minus->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3);

  // Setting Forward Detector acceptance cut
  if((Kaon_Plus->Theta()*TMath::RadToDeg()>5 && Kaon_Plus->Theta()*TMath::RadToDeg()<35)){
    if((Proton->Theta()*TMath::RadToDeg()>5 && Proton->Theta()*TMath::RadToDeg()<35)){
      if((Pi_Minus->Theta()*TMath::RadToDeg()>5 && Pi_Minus->Theta()*TMath::RadToDeg()<35)){
        h_Beam_Momentum2->Fill(Beam.Rho());
        h_Sigma_Momentum_Constrained->Fill(Sigma_Minus->Rho());

        // Filling angular distribution histograms (x,y,weights applied)
        //Vertex 1
        h_cKaon_Plus->Fill(Kaon_Plus->Rho(),Kaon_Plus->Theta()*TMath::RadToDeg(),Phasespace_Weight_1);
        h_cSigma_Minus->Fill(Sigma_Minus->Rho(),Sigma_Minus->Theta()*TMath::RadToDeg(),Phasespace_Weight_1);
        //Vertex 2
        h_cLambda->Fill(Lambda->Rho(),Lambda->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2);
        h_cNeutron->Fill(Neutron->Rho(),Neutron->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2);
        //Vertex 3
        h_cProton->Fill(Proton->Rho(),Proton->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3);
        h_cPi_Minus->Fill(Pi_Minus->Rho(),Pi_Minus->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3);
      }
    }
  }
  // Path Lengths
  h_l_xy->Fill(Phi_Sigma,l_xy);
  h_l_z->Fill(Theta_Sigma*TMath::RadToDeg(),l_z);
  h_Path_Length->Fill(Theta_Sigma*TMath::RadToDeg(),Path_Length);
  h_Path_Length_Decay->Fill(Theta_Sigma*TMath::RadToDeg(),Path_Length_Decay);
}
////////////////////////////////////////////////////////// End of Loop  //////////////////////////////////////////////////////////
h_Beam_Momentum1->Rebin(100);
h_Beam_Momentum2->Rebin(100);
h_Beam_Momentum3->Rebin(100);
h_Beam_Momentum3->Divide (h_Beam_Momentum2,h_Beam_Momentum1);
h_Sigma_Acceptance->Divide(h_Sigma_Momentum_Constrained,h_Sigma_Momentum);
Double_t Scale1 = 1.7577037*100000000/h_Sigma_Momentum1->Integral();
h_Sigma_Momentum1->Scale(Scale1);

for (Int_t i1=1;i1<12001;i1++){
  h_Sigma_Momentum2->SetBinContent(i1,20000);
}

h_Sigma_Scattered->Multiply(h_Sigma_Momentum1,h_Sigma_Momentum2);
h_Sigma_Scattered->Scale((h_Filled_Path_Length->GetMean())*0.16*3/1000000);
h_Sigma_Detected->Multiply(h_Sigma_Scattered,h_Sigma_Acceptance);

cout << "Photon Flux:" << "       " << h_Beam_Flux->Integral() << "  Photons per second" << endl;
cout << "Sigma Flux:" << "        " << h_Sigma_Acceptance->Integral() << "  Sigma particles per second" << endl;
cout << "Sigma Flux:" << "        " << h_Sigma_Acceptance->Integral()*100*24*60*60 << "  Sigma particles per 100 days" << endl;
cout << "Sigma Rescattered:" << " " << h_Sigma_Scattered->Integral() << "  Sigma particles rescattered" << endl;
cout << "Sigma detected:" << "    " << h_Sigma_Detected->Integral() << "  Sigma particles detected in 100 days" << endl;

//////////////////// Rebinning Histograms to make them look more streamlined ////////////////////
h_Sigma_Scattered->Rebin(100);
h_Sigma_Momentum1->Rebin(100);
h_Sigma_Momentum_Constrained->Rebin(100);
h_Sigma_Momentum->Rebin(100);
h_Sigma_Acceptance->Rebin(100);
h_Sigma_Detected->Rebin(100);

//Writes histograms to file
fileOutput1.Write();

////////////////////////////////////////////////////////// Draw Plots //////////////////////////////////////////////////////////
TCanvas* c1 = new TCanvas("c1","Angular Dependence",5,5,1000,1000);
c1->Divide(3,2);
TCanvas* c2 = new TCanvas("c2","Constrained Angular Dependence",5,5,1000,1000);
c2->Divide(3,2);
TCanvas* c3 = new TCanvas("c3","Path Length",5,5,1000,1000);
c3->Divide(3,2);
TCanvas *c4 = new TCanvas("c4","Momentum",5,5,1000,1000);
c4->Divide(2,1);
TCanvas* c5 = new TCanvas("c5","Number of Scattering Events",5,5,1000,1000);
c5->Divide(1,1);
TCanvas* c6 = new TCanvas("c6","Other Plots",5,5,1000,1000);
c6->Divide(3,1);
////////////////// Angular Dependence //////////////////
c1->cd(1);
h_Kaon_Plus->Draw("colz");
c1->cd(2);
h_Sigma_Minus->Draw("colz");
c1->cd(3);
h_Lambda->Draw("colz");
c1->cd(4);
h_Neutron->Draw("colz");
c1->cd(5);
h_Proton->Draw("colz");
c1->cd(6);
h_Pi_Minus->Draw("colz");
////////////////// Constrained Angular Dependence //////////////////
c2->cd(1);
h_cKaon_Plus->Draw("colz");
c2->cd(2);
h_cSigma_Minus->Draw("colz");
c2->cd(3);
h_cLambda->Draw("colz");
c2->cd(4);
h_cNeutron->Draw("colz");
c2->cd(5);
h_cProton->Draw("colz");
c2->cd(6);
h_cPi_Minus->Draw("colz");
////////////////// Path Length //////////////////
c3->cd(1);
h_l_xy->Draw("colz");
c3->cd(2);
h_l_z->Draw("colz");
c3->cd(3);
h_Path_Length->Draw("colz");
c3->cd(4);
h_Path_Length_Decay->Draw("colz");
c3->cd(5);
Decay_Plot->Draw();
c3->cd(6);
h_Filled_Path_Length->Draw();
////////////////// Momentum //////////////////
c4->cd(1);
h_Sigma_Momentum->Draw();
c4->cd(2);
h_Beam_Momentum6->Draw("hist");
////////////////// Scattering Events //////////////////
c5->cd(1);
h_Sigma_Scattered->Draw();
////////////////// Others //////////////////
c6->cd(1);
h_Beam_Flux->Draw("colz");
c6->cd(2);
h_Sigma_Acceptance->Draw("colz");
c6->cd(3);
h_Sigma_Detected->Draw("colz");
}