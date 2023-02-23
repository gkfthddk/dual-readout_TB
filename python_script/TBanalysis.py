import ROOT
import pydrcTB
import ctypes

chain = ROOT.TChain("events")
chain.Add("ntuple_Run_678_Wave_FILE_0.root")
# chain.Add("ntuple_Run_642_Fast_FILE_0.root")
# chain.Add("ntuple_Run_643_Fast_FILE_0.root")

utility = pydrcTB.TButility()
utility.loading("mapping_data_MCPPMT_positiveSignal_v3.csv")

channelsize = 32

output = ROOT.TFile("output.root","RECREATE")

# histoPS = ROOT.TH1D("histoPS","PS adc",1000,-10000,90000)
# # histoPSnoise = ROOT.TH1D("histoPSnoise","first 50 bin ADC - last 50 bin ADC",800,-5000,3000)
# histoSiPMadc_S = ROOT.TH1D("histoSiPMadc_S","SiPM ADC sum S",1000,-10000,190000)
# histoSiPMadc_C = ROOT.TH1D("histoSiPMadc_C","SiPM ADC sum C",1000,-10000,190000)
# histoSiPMnoise = ROOT.TH1D("histoSiPMnoise","first 50 bin ADC - last 50 bin ADC",1000,-10000,19000)
# histoMuon = ROOT.TH1D("histoMuon","Muon ADC",1000,-5000,45000)
histoSiPM_C = ROOT.TH1D("histoSiPM_C","individual SiPM ADC C",1000,-5000,95000)
histoSiPM_S = ROOT.TH1D("histoSiPM_S","individual SiPM ADC S",1000,-5000,95000)

def checkAux(anevt,atype):
    for imid in range(anevt.size()):
        for ich in range(channelsize):
            cidboost = pydrcTB.TBcid( imid+1, ich+1 )
            adet = utility.find(cidboost)

            if (not adet.isNull()) and (not adet.isModule()) and adet.detType()==atype:
                return imid+1, ich+1

    raise RuntimeError('checkPS: No PS detector found!')

PSmid = -1
PSch = -1
MuonMid = -1
MuonCh = -1

for ievt in range(chain.GetEntries()):
    chain.GetEntry(ievt)
    anevt = getattr(chain,"TBevt")

    if ievt>1000:
        break

    if ievt==0:
        PSmid, PSch = checkAux(anevt,25)
        print('(%d, %d)' % (PSmid,PSch))
        MuonMid, MuonCh = checkAux(anevt,27)
        print('(%d, %d)' % (MuonMid,MuonCh))

    if ievt % 10000 == 0:
        print("Processing %d th event" % (ievt))

    cidPS = ROOT.TBcid( PSmid, PSch )
    dataPS = anevt.data( cidPS )

    # adc_PS = dataPS.emulfastADC(40,120)
    # histoPSnoise.Fill(noiseLevel(awave))
    # histoPS.Fill( adata.adc() )

    # cidMuon = ROOT.TBcid( MuonMid, MuonCh )
    # adc_mu = anevt.data(cidMuon).adc()

    # adc_S = 0.
    # adc_C = 0.

    # if adc_PS < 25000.:
    #     continue
    #
    # if adc_mu > 3000:
    #     continue

    for imid in range(anevt.size()):
        for ich in range(channelsize):
            cidboost = pydrcTB.TBcid( imid+1, ich+1 )
            cidroot = ROOT.TBcid( imid+1, ich+1 )
            adet = utility.find(cidboost)
            adata = anevt.data(cidroot)
            # awave = adata.waveform()

            if (not adet.isNull()) and adet.isModule() and adet.isSiPM() and adet.isCeren():
                # adc_C += adata.emulfastADC(40,120)
                # histoSiPM_C.Fill(adata.adc())

                if adata.emulfastADC(40,120) > 2000.:
                    histo = ROOT.TH1D("C_histo"+"_evt"+str(ievt)+"_mid"+str(imid+1)+"_pl"+str(adet.plate())+"_col"+str(adet.column()),"C_waveform"+"_evt"+str(ievt)+"_mid"+str(imid+1)+"_pl"+str(adet.plate())+"_col"+str(adet.column()),1000,0,1000)

                    awave = adata.waveform()

                    max = 0

                    for bin in range(1,1001):
                        histo.SetBinContent(bin,awave[bin])

                        if max < awave[bin]:
                            max = awave[bin]

                    for bin in range(1,1001):
                        val = max - histo.GetBinContent(bin)
                        histo.SetBinContent(bin,val)

                    output.WriteTObject(histo)

                    ROOT.TVirtualFFT.SetTransform(0)

                    histo_fft = 0
                    histo_fft = histo.FFT(histo_fft,"MAG M")

                    histo_fft_fi = ROOT.TH1D("C_fft"+"_evt"+str(ievt)+"_mid"+str(imid+1)+"_pl"+str(adet.plate())+"_col"+str(adet.column()),"C_fft"+"_evt"+str(ievt)+"_mid"+str(imid+1)+"_pl"+str(adet.plate())+"_col"+str(adet.column()),histo_fft.GetNbinsX(),0,histo_fft.GetNbinsX())

                    for bin in range(1,1001):
                        histo_fft_fi.SetBinContent(bin,histo_fft.GetBinContent(bin))

                    del histo_fft

                    output.WriteTObject(histo_fft_fi)

                    nbins = int(histo_fft_fi.GetNbinsX())
                    nbins_new = int(histo_fft_fi.GetNbinsX()/2+1)

                    fft = ROOT.TVirtualFFT.GetCurrentTransform()
                    re = []
                    im = []

                    for bin in range(0,histo.GetNbinsX()):
                        re_i = ctypes.c_double(0.)
                        im_i = ctypes.c_double(0.)
                        fft.GetPointComplex(bin,re_i,im_i)
                        re.append(re_i)
                        im.append(im_i)

                    for bin in range(1,5):
                        re[bin] = 0.
                        im[bin] = 0.
                    #
                    for bin in range(490,510):
                        re[bin] = 0.
                        im[bin] = 0.

                    nbin_cpp = ctypes.c_int(nbins)

                    fft_own = ROOT.TVirtualFFT.FFT(1,nbin_cpp,"C2R M")

                    for bin in range(0,nbins_new):
                        fft_own.SetPoint(bin,re[bin],im[bin])

                    fft_own.Transform()

                    ans = 0
                    ans = ROOT.TH1.TransformHisto(fft_own,ans,"MAG M")

                    histo_rev_fi = ROOT.TH1D("C_rev"+"_evt"+str(ievt)+"_mid"+str(imid+1)+"_pl"+str(adet.plate())+"_col"+str(adet.column()),"C_rev"+"_evt"+str(ievt)+"_mid"+str(imid+1)+"_pl"+str(adet.plate())+"_col"+str(adet.column()),1000,0,1000)

                    for bin in range(0,1000):
                        histo_rev_fi.SetBinContent(bin,ans.GetBinContent(bin))

                    del ans

                    output.WriteTObject(histo_rev_fi)

                # adc_C += anevt.data( cidroot ).adc()
            elif (not adet.isNull()) and adet.isModule() and adet.isSiPM() and (not adet.isCeren()):
                if adata.emulfastADC(40,120) > 2000.:
                    histo = ROOT.TH1D("S_histo"+"_evt"+str(ievt)+"_mid"+str(imid+1)+"_pl"+str(adet.plate())+"_col"+str(adet.column()),"S_waveform"+"_evt"+str(ievt)+"_mid"+str(imid+1)+"_pl"+str(adet.plate())+"_col"+str(adet.column()),1000,0,1000)

                    awave = adata.waveform()

                    max = 0

                    for bin in range(1,1001):
                        histo.SetBinContent(bin,awave[bin])

                        if max < awave[bin]:
                            max = awave[bin]

                    for bin in range(1,1001):
                        val = max - histo.GetBinContent(bin)
                        histo.SetBinContent(bin,val)

                    output.WriteTObject(histo)

                    ROOT.TVirtualFFT.SetTransform(0)

                    histo_fft = 0
                    histo_fft = histo.FFT(histo_fft,"MAG M")

                    histo_fft_fi = ROOT.TH1D("S_fft"+"_evt"+str(ievt)+"_mid"+str(imid+1)+"_pl"+str(adet.plate())+"_col"+str(adet.column()),"S_fft"+"_evt"+str(ievt)+"_mid"+str(imid+1)+"_pl"+str(adet.plate())+"_col"+str(adet.column()),histo_fft.GetNbinsX(),0,histo_fft.GetNbinsX())

                    for bin in range(1,1001):
                        histo_fft_fi.SetBinContent(bin,histo_fft.GetBinContent(bin))

                    del histo_fft

                    output.WriteTObject(histo_fft_fi)

                    nbins = int(histo_fft_fi.GetNbinsX())
                    nbins_new = int(histo_fft_fi.GetNbinsX()/2+1)

                    fft = ROOT.TVirtualFFT.GetCurrentTransform()
                    re = []
                    im = []

                    for bin in range(0,histo.GetNbinsX()):
                        re_i = ctypes.c_double(0.)
                        im_i = ctypes.c_double(0.)
                        fft.GetPointComplex(bin,re_i,im_i)
                        re.append(re_i)
                        im.append(im_i)

                    for bin in range(1,5):
                        re[bin] = 0.
                        im[bin] = 0.
                    #
                    for bin in range(490,510):
                        re[bin] = 0.
                        im[bin] = 0.

                    nbin_cpp = ctypes.c_int(nbins)

                    fft_own = ROOT.TVirtualFFT.FFT(1,nbin_cpp,"C2R M")

                    for bin in range(0,nbins_new):
                        fft_own.SetPoint(bin,re[bin],im[bin])

                    fft_own.Transform()

                    ans = 0
                    ans = ROOT.TH1.TransformHisto(fft_own,ans,"MAG M")

                    histo_rev_fi = ROOT.TH1D("S_rev"+"_evt"+str(ievt)+"_mid"+str(imid+1)+"_pl"+str(adet.plate())+"_col"+str(adet.column()),"S_rev"+"_evt"+str(ievt)+"_mid"+str(imid+1)+"_pl"+str(adet.plate())+"_col"+str(adet.column()),1000,0,1000)

                    for bin in range(0,1000):
                        histo_rev_fi.SetBinContent(bin,ans.GetBinContent(bin))

                    del ans

                    output.WriteTObject(histo_rev_fi)

            #     histoSiPM_S.Fill(adata.adc())
            #     adc_S += adata.emulfastADC(40,120)
            #     histoSiPMnoise.Fill(noise)
                # adc_S += anevt.data( cidroot ).adc()
            # elif (not adet.isNull()) and (not adet.isModule()) and adet.detType()==27:
                # adc_mu = anevt.data( cidroot ).adc()

    # histoSiPMadc_S.Fill(adc_S)
    # histoSiPMadc_C.Fill(adc_C)
    # histoMuon.Fill(adc_mu)
    # histoPS.Fill(adc_PS)

# output.WriteTObject(histoPS)
# output.WriteTObject(histoMuon)
# output.WriteTObject(histoPSnoise)
# output.WriteTObject(histoSiPMadc_S)
# output.WriteTObject(histoSiPMadc_C)
# output.WriteTObject(histoSiPMnoise)
# output.WriteTObject(histoSiPM_C)
# output.WriteTObject(histoSiPM_S)
output.Close()
