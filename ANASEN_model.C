#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"

#include "TMath.h"

void ANASEN_model() {

  // Create ROOT manager and master volume

  TGeoManager *geom = new TGeoManager("Detector", "3D Detector");

  //--- define some materials
  TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
  TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
  //--- define some media
  TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
  TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);

  //--- make the top container volume 
  Double_t worldx = 100.; //cm
  Double_t worldy = 100.; //cm
  Double_t worldz = 100.; //cm
  TGeoVolume *worldBox = geom->MakeBox("ROOT", Vacuum, worldx, worldy, worldz);
  geom->SetTopVolume(worldBox);


  TGeoVolume *axisX = geom->MakeTube("axisX", Al, 0, 0.1, 5.);
  axisX->SetLineColor(1);
  worldBox->AddNode(axisX, 1, new TGeoCombiTrans(5, 0, 0., new TGeoRotation("rotA", 90., 90., 0.)));

  TGeoVolume *axisY = geom->MakeTube("axisY", Al, 0, 0.1, 5.);
  axisY->SetLineColor(1);
  worldBox->AddNode(axisY, 1, new TGeoCombiTrans(0, 5, 0., new TGeoRotation("rotB", 0., 90., 0.)));

  TGeoVolume *axisZ = geom->MakeTube("axisZ", Al, 0, 0.1, 5.);
  axisZ->SetLineColor(1);
  worldBox->AddNode(axisZ, 1, new TGeoTranslation(0, 0,  5));

  const int nWire = 24;
  const int zLen = 15;
  TGeoVolume *pcA = geom->MakeTube("tub1", Al, 0, 0.01, zLen);
  pcA->SetLineColor(4);
  for( int i = 0; i < nWire; i++){
    worldBox->AddNode(pcA, i+1, new TGeoCombiTrans( 3.8 * TMath::Cos( TMath::TwoPi() / nWire *i), 
                                                    3.8 * TMath::Sin( TMath::TwoPi() / nWire *i), 
                                                    0, 
                                                    new TGeoRotation("rot1", 360/nWire * (i), 20., 0.)));
  }


  TGeoVolume *pcC = geom->MakeTube("tub2", Al, 0, 0.01, zLen);
  pcC->SetLineColor(6);
  for( int i = 0; i < nWire; i++){
    worldBox->AddNode(pcC, i+1, new TGeoCombiTrans( 4.3 * TMath::Cos( TMath::TwoPi() / nWire *i), 
                                                    4.3 * TMath::Sin( TMath::TwoPi() / nWire *i), 
                                                    0, 
                                                    new TGeoRotation("rot1", 360/nWire * (i), -20., 0.)));
  }

  const int nSX3 = 12;
  const int sx3Radius = 8.8;
  TGeoVolume * sx3 = geom->MakeBox("box", Al, 0.1, 4.0/2, 7.5/2);
  sx3->SetLineColor(2);
  for( int i = 0; i < nSX3; i++){
    worldBox->AddNode(sx3, 2*i+1., new TGeoCombiTrans( sx3Radius * TMath::Cos( TMath::TwoPi() / nSX3 *i), 
                                                       sx3Radius * TMath::Sin( TMath::TwoPi() / nSX3 *i), 
                                                     8.0/2, 
                                                     new TGeoRotation("rot1", 360/nSX3 * (i), 0., 0.)));

    worldBox->AddNode(sx3, 2*i+2., new TGeoCombiTrans( sx3Radius* TMath::Cos( TMath::TwoPi() / nSX3 *i), 
                                                     sx3Radius * TMath::Sin( TMath::TwoPi() / nSX3 *i), 
                                                     -8.0/2, 
                                                     new TGeoRotation("rot1", 360/nSX3 * (i), 0., 0.)));
  }

  geom->CloseGeometry();

  geom->SetVisLevel(4);
  worldBox->Draw("ogle");
}
