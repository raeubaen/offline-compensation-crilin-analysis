import argparse
import ROOT
import array

def setStyle():
    ROOT.gStyle.SetTitleFont(132,"xyz");
    ROOT.gStyle.SetTitleFont(132," ");
    ROOT.gStyle.SetTitleSize(0.06,"xyz");
    ROOT.gStyle.SetTitleSize(0.06," ");
    ROOT.gStyle.SetLabelFont(132,"xyz");
    ROOT.gStyle.SetLabelSize(0.05,"xyz");
    ROOT.gStyle.SetTextFont(132);
    ROOT.gStyle.SetTextSize(0.08);
    ROOT.gStyle.SetStatFont(132);

def makeCut(cuts,verbose=False):
    cuts = ["("+c+")" for c in cuts]
    cut="*".join(cuts)
    if verbose:
        print("Selection to be applied = ",cut)
    return cut

def gaussFit(h, name, title, xmin=-1, xmax=-1):

    # --- Create RooFit variable ---
    x = ROOT.RooRealVar("x", "E/E_{True}",
                        h.GetXaxis().GetXmin(),
                        h.GetXaxis().GetXmax())

    # --- Define fit subrange (optional) ---
    useRange = False
    if xmin >= 0 and xmax >= 0:
        x.setRange("fitRange", xmin, xmax)
        useRange = True

    # --- Import histogram ---
    data = ROOT.RooDataHist("data", "data",
                            ROOT.RooArgList(x), h)

    # --------------------------------
    # Gaussian parameters
    # --------------------------------
    mean  = ROOT.RooRealVar("mean", "Gaussian mean",
                            h.GetMean(),
                            h.GetMean()-1,
                            h.GetMean()+1)

    sigma = ROOT.RooRealVar("sigma", "Gaussian sigma",
                            h.GetRMS(),
                            0.1*h.GetRMS(),
                            5*h.GetRMS())

    # --------------------------------
    # Gaussian PDF
    # --------------------------------
    gauss = ROOT.RooGaussian("gauss", "Gaussian",
                             x, mean, sigma)

    # --------------------------------
    # Extended yield
    # --------------------------------
    nsig = ROOT.RooRealVar("nsig", "signal yield",
                           h.Integral(),
                           0.0,
                           10.0*h.Integral())

    model = ROOT.RooAddPdf("model", "extended Gaussian model",
                           ROOT.RooArgList(gauss),
                           ROOT.RooArgList(nsig))

    # ---------------------------
    # Fit
    # ---------------------------
    fitArgs = [ROOT.RooFit.Extended(True),
               ROOT.RooFit.Save(),
               ROOT.RooFit.PrintLevel(-1)]

    if useRange:
        fitArgs.insert(0, ROOT.RooFit.Range("fitRange"))

    result = model.fitTo(data, *fitArgs)

    # ---------------------------
    # Plot
    # ---------------------------
    frame = x.frame(ROOT.RooFit.Title(title))
    data.plotOn(frame)
    model.plotOn(frame)

    c = ROOT.TCanvas("c", "Gaussian Fit", 800, 600)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)

    frame.GetYaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.Draw()

    # Info box
    pt = ROOT.TPaveText(0.60, 0.70, 0.88, 0.88, "NDC")
    pt.SetFillColor(0)
    pt.SetTextFont(42)
    pt.SetBorderSize(0)
    pt.SetTextSize(0.05)

    pt.AddText(f"#mu = {mean.getVal():.3g} #pm {mean.getError():.3g}")
    pt.AddText(f"#sigma = {sigma.getVal():.3g} #pm {sigma.getError():.3g}")

    pt.Draw()

    for ext in ["png", "pdf"]:
        c.SaveAs(f"{name}.{ext}")

    print("Fit results:")
    result.Print()

    return {"mean": (mean.getVal(), mean.getError()),
            "sigma": (sigma.getVal(), sigma.getError())}


def cbFit(h,name,title,xmin=-1,xmax=-1):
    # --- Create RooFit variables ---
    x = ROOT.RooRealVar("x", "E/E_{True}", h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())

    # --- Define the fit subrange ---
    if xmin>=0 and xmax>=0:
        x.setRange("fitRange", xmin, xmax)

    # --- Import histogram into RooDataHist ---
    data = ROOT.RooDataHist("data", "data", ROOT.RooArgList(x), h)

    # --------------------------------
    # Double Crystal Ball parameters
    # --------------------------------
    mean  = ROOT.RooRealVar("mean",  "DCB mean",
                            h.GetMean(),
                            h.GetMean()-1,
                            h.GetMean()+1)

    sigma = ROOT.RooRealVar("sigma", "DCB sigma",
                            h.GetRMS(),
                            0.1*h.GetRMS(),
                            5*h.GetRMS())

    alphaL = ROOT.RooRealVar("alphaL", "alphaL", 1.5, 0.1, 5.0)
    nL     = ROOT.RooRealVar("nL",     "nL",     3.0, 0.5, 20.0)

    alphaR = ROOT.RooRealVar("alphaR", "alphaR", 1.5, 0.1, 5.0)
    nR     = ROOT.RooRealVar("nR",     "nR",     3.0, 0.5, 20.0)

    # --------------------------------
    # Double Crystal Ball PDF
    # --------------------------------
    dcb = ROOT.RooCrystalBall(
        "dcb", "Double Crystal Ball",
        x,
        mean,
        sigma,
        alphaL, nL,
        alphaR, nR
     )

    # --------------------------------
    # Extended yield
    # --------------------------------
    nsig = ROOT.RooRealVar("nsig", "signal yield",
                           h.Integral(),
                           0.0,
                           10.0*h.Integral())

    # --------------------------------
    # Final model = ONLY DCB
    # --------------------------------
    model = ROOT.RooAddPdf("model", "extended DCB model",
                           ROOT.RooArgList(dcb),
                           ROOT.RooArgList(nsig))

    # ---------------------------
    # Fit (extended likelihood)
    # ---------------------------
    result = model.fitTo(data,
                          ROOT.RooFit.Range("fitRange"),
                          ROOT.RooFit.Extended(True),
                          ROOT.RooFit.Save(),
                          ROOT.RooFit.PrintLevel(-1))


    # ---------------------------
    # Plot
    # ---------------------------
    frame = x.frame(ROOT.RooFit.Title(f"{title}"))
    data.plotOn(frame)
    model.plotOn(frame)
    
    # Individual components
    model.plotOn(frame, ROOT.RooFit.Components("gauss"),
                 ROOT.RooFit.LineColor(ROOT.kBlue),
                 ROOT.RooFit.LineStyle(ROOT.kDashed))
    
    model.plotOn(frame, ROOT.RooFit.Components("cb"),
                 ROOT.RooFit.LineColor(ROOT.kRed),
                 ROOT.RooFit.LineStyle(ROOT.kDashed))
    
    c = ROOT.TCanvas("c", "Extended Fit", 800, 600)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetTitleOffset(1.1)

    frame.Draw()
    
    # Add info box with mean and sigma
    pt = ROOT.TPaveText(0.60, 0.65, 0.88, 0.88, "NDC")
    pt.SetFillColor(0)
    pt.SetTextFont(42)
    pt.SetBorderSize(0)
    pt.SetTextSize(0.05)

    pt.AddText(f"m_{{core}} = {mean.getVal():.3g} #pm {mean.getError():.3g}")
    pt.AddText(f"#sigma_{{core}} = {sigma.getVal():.3g} #pm {sigma.getError():.3g}")
    
    pt.Draw()
    
    for ext in ["png","pdf","root"]:
        c.SaveAs(f"{name}.{ext}")

    print("Fit results:")
    result.Print()
    print("SIGMA: ", sigma.getVal())
    return {"mean":(mean.getVal(),mean.getError()),
            "sigma":(sigma.getVal(),sigma.getError())}

def plotSingleResolution(tree,name,title,selection=[],verbose=False,nbins=200,just_gauss=False,time=False):
    sel=makeCut(selection,verbose)

    canvas = ROOT.TCanvas("c1", "resolutions", 800, 600)

    tree.Draw("Sum$(Hit_NCherenkov / 24.3)/(PrimaryEnergy) - 1>>resotemp(1000, 0, 1)", "", "goff")

    resotemp = ROOT.gDirectory.Get("resotemp")
    m = resotemp.GetMean()
    s = resotemp.GetRMS()
    reso = ROOT.TH1F("reso","resolution",nbins,m-5*s,m+5*s)
    tree.Draw("Sum$(Hit_NCherenkov / 24.3)/(PrimaryEnergy) - 1>> reso",sel, "goff")
    if not just_gauss: results = cbFit(reso,name,title)
    else: results = gaussFit(reso,name,title)
    return results

def plotDifferentialResolution(tree,selection=[],verbose=False,time=False):
    Ebins=list(range(0, 100, 5))
    x,ex,b,eb,s,es = [],[],[],[],[],[]
    print ("Energy bins to be analysed: ",Ebins)
    for ie in range(len(Ebins)-1):
        print ("ie = ",ie)
        addcut = f"PrimaryEnergy[0]>{Ebins[ie]} && PrimaryEnergy[0]<={Ebins[ie+1]}"
        print ("Processing bin: ",addcut)
        fullsel = selection + [addcut]

        name = f"resolution_E{Ebins[ie]}To{Ebins[ie+1]}"
        title = f"{Ebins[ie]} GeV < E < {Ebins[ie+1]} GeV"
        results = plotSingleResolution(tree,name,title,fullsel, nbins=100, just_gauss=False, time=time)
        b.append(results["mean"][0])
        eb.append(results["mean"][1])
        s.append(results["sigma"][0])
        es.append(results["sigma"][1])
        x.append((Ebins[ie]+Ebins[ie+1])/2.)
        ex.append(0)

    print("bias, errorbias, sigma, errorsigma, x, errorx", b,eb,s,es,x,ex)

    # Convert Python lists to C-style arrays
    x_arr = array.array('d', x)
    ex_arr = array.array('d', ex)
    b_arr = array.array('d', b)
    eb_arr = array.array('d', eb)
    s_arr = array.array('d', s)
    es_arr = array.array('d', es)

    gbias = ROOT.TGraphErrors(len(x), x_arr, b_arr, ex_arr, eb_arr)
    gsigma = ROOT.TGraphErrors(len(x), x_arr, s_arr, ex_arr, es_arr)

    # Optional: style
    gbias.SetTitle("")
    gbias.SetMarkerStyle(20)
    gbias.GetXaxis().SetTitle("Amplitude [GeV]")
    gbias.GetYaxis().SetTitle("Mean of E/E_{True} - 1")
    
    gsigma.SetTitle("")
    gsigma.SetMarkerStyle(20)
    gsigma.GetXaxis().SetTitle("Amplitude [GeV]")
    gsigma.GetYaxis().SetTitle("#sigma of [E/E_{True} - 1]")

    # Draw
    c = ROOT.TCanvas("c", "TGraphErrors from arrays", 800, 600)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    gbias.GetYaxis().SetTitleOffset(1.1)
    gbias.GetXaxis().SetTitleOffset(1.1)
    gsigma.GetYaxis().SetTitleOffset(1.1)
    gsigma.GetXaxis().SetTitleOffset(1.1)
    gsigma.GetYaxis().SetRangeUser(0,0.20)

    gbias.Draw("APC")
    for ext in ["png","pdf","root"]:
        c.SaveAs(f"bias_vs_E.{ext}")
    gsigma.Draw("APC")
    for ext in ["png","pdf","root"]:
        c.SaveAs(f"sigma_vs_E.{ext}")
    
def main():
    parser = argparse.ArgumentParser(
        description="Script to run simple analysis from a reconstructed multifit TTree"
    )

    # Required positional argument
    parser.add_argument("input_file", help="Path to the input file")

    # Optional flags
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Enable verbose mode")

    # Option that can appear multiple times
    parser.add_argument("-c", "--cut",
                        action="append",
                        default=["1"],
                        help="cut to be added (can be given multiple times with list syntax)")

    parser.add_argument("-a", "--analysis",
                        action="append",
                        default=[],
                        help="Analysis that should be run (can be given multiple times with list syntax): single_resolution or differential_resolution")


    parser.add_argument("-o", "--output",
                        type=str,
                        help="Output file name")

    args = parser.parse_args()

    tfile = ROOT.TFile.Open(args.input_file)
    tree = tfile.Get("events")

    tree.SetAlias("rng", "sin(2*pi*rndm)*sqrt(-2*log(rndm))")

    setStyle()

    print (" ==== Analyses to be run: ====\n", args.analysis)

    if "single_resolution" in args.analysis:
        plotSingleResolution(tree,"resolution"," AND ".join(args.cut),args.cut,args.verbose)

    if "differential_resolution" in args.analysis:
        plotDifferentialResolution(tree,args.cut,args.verbose)

if __name__ == "__main__":
    main()
