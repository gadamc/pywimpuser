#!/usr/bin/env python

from pyWIMP.DMModels.base_model import BaseVariables
import pyWIMP.DMModels.wimp_model as wimp_model
import ROOT
import re

#pair of mass / sigmas to plot
mass_sigmas = ([7, 1.0e-4], [9, 8.0e-5], [10, 2.0e-5], [30, 1.0e-5])

minEnergy = 0
energyBinSize = 50
nuclRecoilScale = True
if nuclRecoilScale: maxEnergy = 20
else: maxEnergy = 5

#nuclear recoil keV scale example
if nuclRecoilScale: energyRanges = ([3.0, 5.0], [5.0, 8.0], [5.0, 7.3], [3.0, 8.0], [5.0, 11.9], [5.0, 15.0])
#electron equivalent keV scale example
else: energyRanges = ([0.5, 0.9], [0.5, 3.0], [0.9, 3.0], [1.1, 2.62])

startTime = 0.0 
totalTime = 2.0 #years

def makeHist(lowE, highE, sigma, binscale):
  hist1d = hist.ProjectionY("counts_%sto%skeV"%(str(lowE), str(highE)), int(lowE*binscale), int(highE*binscale))
  hist1d.Scale(1.0/(binscale * (highE - lowE)) )
  hist1d.SetTitle("Event Rate in %s to %s keV range at #sigma = %s pb" % (str(lowE), str(highE), sigma))
  hist1d.SetYTitle("#frac{dR}{dE} [counts/keV/kg/day] ")
  hist1d.SetXTitle("Time since Jan 1 [years]")
  hist1d.Write()
  
for mspair in mass_sigmas:
  mass_of_wimp = mspair[0]
  sigma_nucleon = mspair[1]
  print mass_of_wimp, sigma_nucleon

  basevars = BaseVariables(startTime, totalTime, minEnergy, maxEnergy, True)  #could add offset, but defaults to Jan 1
  time = basevars.get_time()
  time.setConstant()
  time.setVal(0)
  
  wm = wimp_model.WIMPModel(basevars, mass_of_wimp,  nucl_recoil=nuclRecoilScale)
  model = wm.get_WIMP_model()
  expected_events = model.expectedEvents(ROOT.RooArgSet(basevars.get_energy()))
  hist = model.createHistogram("hist", basevars.get_energy(), \
                          ROOT.RooFit.Binning( int(maxEnergy*energyBinSize) ),\
                          ROOT.RooFit.YVar(basevars.get_time(), \
                          ROOT.RooFit.Binning(int(totalTime*100))))
  hist.SetLineColor(4)

  A = wm.atomic_mass_of_target
  mass_of_target = wm.mass_of_target.getVal(ROOT.RooArgSet(A))
  mu_A = mass_of_wimp * mass_of_target/ (mass_of_wimp + mass_of_target)
  mu_one = 0.925
  sigma_nucleus =  sigma_nucleon * A.getVal()**2 * mu_A**2 / mu_one**2 

  hist.Scale(sigma_nucleus)
  hist.Scale(1./365.25) #scale years to day
  hist.GetXaxis().CenterTitle()
  if nuclRecoilScale:
    energyUnit = 'keVnr'
  else: 
    energyUnit = 'keVee'
    
  hist.GetXaxis().SetTitle("Energy [%s]" % (energyUnit,))
  hist.GetXaxis().SetTitleOffset(1.5)
  hist.GetYaxis().CenterTitle()
  hist.GetYaxis().SetTitle("Time [years]")
  hist.GetYaxis().SetTitleOffset(1.5)
  hist.GetYaxis().SetNdivisions(507)
  hist.GetZaxis().CenterTitle()
  hist.GetZaxis().SetTitle("Amplitude [a.u.]")
  hist.GetZaxis().SetTitle("#frac{dR}{dE} [counts/keV/kg/day]")
  hist.GetZaxis().SetTitleOffset(1.1)
  hist.GetZaxis().SetNdivisions(509)
  ysize = hist.GetYaxis().GetBinWidth(1)
  xsize = hist.GetXaxis().GetBinWidth(1)
  time = basevars.get_time()
  time_scale = time.getMax() - time.getMin()
  print ysize, xsize
  hist.Scale(expected_events/(xsize*ysize))
  c1 = ROOT.TCanvas()
  c1.SetLogz()
  c1.SetTheta(22.1822);
  c1.SetPhi(-24.31034);
  hist.Draw("SURF4FB")
  hist.GetZaxis().SetRangeUser(0.05, 1000)
  c1.Update()
  rootFile = '%dGeV_%s.root' % (mass_of_wimp, energyUnit)
  print rootFile
  f = ROOT.TFile(rootFile,'recreate')
  hist.SetName("3dhist")
  hist.Write()
  
  for e_range in energyRanges:
    makeHist(e_range[0],e_range[1],str(sigma_nucleon),energyBinSize)

  c1.Print(rootFile.split('.root')[0]+".eps")
  f.Close()
  
