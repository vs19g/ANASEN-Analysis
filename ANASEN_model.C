#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"

#include "TMath.h"

void ANASEN_model(int anodeID1 = -1, int anodeID2 = -1, int cathodeID1 = -1, int cathodeID2 = -1 ) {

  // Create ROOT manager and master volume
  TGeoManager *geom = new TGeoManager("Detector", "ANASEN");

  //--- define some materials
  TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
  TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
  //--- define some media
  TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
  TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);

  //--- make the top container volume 
  Double_t worldx = 200.; //mm
  Double_t worldy = 200.; //mm
  Double_t worldz = 200.; //mm
  TGeoVolume *worldBox = geom->MakeBox("ROOT", Vacuum, worldx, worldy, worldz);
  geom->SetTopVolume(worldBox);

  //--- making axis
  TGeoVolume *axisX = geom->MakeTube("axisX", Al, 0, 0.1, 5.);
  axisX->SetLineColor(1);
  worldBox->AddNode(axisX, 1, new TGeoCombiTrans(5, 0, 0., new TGeoRotation("rotA", 90., 90., 0.)));

  TGeoVolume *axisY = geom->MakeTube("axisY", Al, 0, 0.1, 5.);
  axisY->SetLineColor(1);
  worldBox->AddNode(axisY, 1, new TGeoCombiTrans(0, 5, 0., new TGeoRotation("rotB", 0., 90., 0.)));

  TGeoVolume *axisZ = geom->MakeTube("axisZ", Al, 0, 0.1, 5.);
  axisZ->SetLineColor(1);
  worldBox->AddNode(axisZ, 1, new TGeoTranslation(0, 0,  5));

  //--- making ANASEN
  const int nWire = 24;
  const int wireShift = 3;
  const int zLen = 300; //mm
  const int radiusA = 38;
  const int radiusC = 43;

  //.......... convert to wire center dimensions
  double dAngle = wireShift * TMath::TwoPi() / nWire;
  double radiusAnew = radiusA * TMath::Cos( dAngle / 2.);
  double wireALength = TMath::Sqrt( zLen*zLen + TMath::Power(2* radiusA * TMath::Sin(dAngle/2),2) );
  double wireATheta = TMath::ATan2( 2* radiusA * TMath::Sin( dAngle / 2.), zLen);

  // printf("    dAngle : %f\n", dAngle);
  // printf(" newRadius : %f\n", radiusAnew);
  // printf("wireLength : %f\n", wireALength);
  // printf("wire Theta : %f\n", wireATheta);

  TGeoVolume *pcA = geom->MakeTube("tub1", Al, 0, 0.01, wireALength/2);
  pcA->SetLineColor(4);  

  for( int i = 0; i < nWire; i++){
    if( i < anodeID1 || i > anodeID2) continue;   
    worldBox->AddNode(pcA, i+1, new TGeoCombiTrans( radiusAnew * TMath::Cos( TMath::TwoPi() / nWire *i + dAngle / 2), 
                                                    radiusAnew * TMath::Sin( TMath::TwoPi() / nWire *i + dAngle / 2), 
                                                    0, 
                                                    new TGeoRotation("rot1", 360/ nWire * (i + wireShift/2.), wireATheta * 180/ TMath::Pi(), 0.)));

  }

  double radiusCnew = radiusC * TMath::Cos( dAngle / 2.);
  double wireCLength = TMath::Sqrt( zLen*zLen + TMath::Power(2* radiusC * TMath::Sin(dAngle/2),2) );
  double wireCTheta = TMath::ATan2( 2* radiusC * TMath::Sin( dAngle / 2.), zLen);

  TGeoVolume *pcC = geom->MakeTube("tub2", Al, 0, 0.01, wireCLength/2);
  pcC->SetLineColor(6);
  for( int i = 0; i < nWire; i++){
    if( i < cathodeID1 || i > cathodeID2) continue;   
    worldBox->AddNode(pcC, i+1, new TGeoCombiTrans( radiusCnew * TMath::Cos( TMath::TwoPi() / nWire *i - dAngle/2), 
                                                    radiusCnew * TMath::Sin( TMath::TwoPi() / nWire *i - dAngle/2), 
                                                    0, 
                                                    new TGeoRotation("rot1", 360/ nWire * (i - wireShift/2.), -wireCTheta * 180/ TMath::Pi(), 0.)));
  }

  const int nSX3 = 12;
  const int sx3Radius = 88;
  const int sx3Width = 40;
  const int sx3Length = 75;
  const int sx3Gap = 5;

  TGeoVolume * sx3 = geom->MakeBox("box", Al, 0.1, sx3Width/2, sx3Length/2);
  sx3->SetLineColor(kGreen+3);
  for( int i = 0; i < nSX3; i++){
    worldBox->AddNode(sx3, 2*i+1., new TGeoCombiTrans( sx3Radius * TMath::Cos( TMath::TwoPi() / nSX3 *i), 
                                                       sx3Radius * TMath::Sin( TMath::TwoPi() / nSX3 *i), 
                                                     sx3Length/2+sx3Gap, 
                                                     new TGeoRotation("rot1", 360/nSX3 * (i), 0., 0.)));

    worldBox->AddNode(sx3, 2*i+2., new TGeoCombiTrans( sx3Radius* TMath::Cos( TMath::TwoPi() / nSX3 *i), 
                                                     sx3Radius * TMath::Sin( TMath::TwoPi() / nSX3 *i), 
                                                     -sx3Length/2-sx3Gap, 
                                                     new TGeoRotation("rot1", 360/nSX3 * (i), 0., 0.)));
  }

  geom->CloseGeometry();

  geom->SetVisLevel(4);
  worldBox->Draw("ogle");
}
