import argparse
import array
import os
import re
import ROOT


def make_cut(cuts, verbose=False):
    cuts = ["("+c+")" for c in cuts]
    cut="*".join(cuts)
    if verbose:
        print("Selection to be applied = ",cut)
    return cut


def reco_expr(kind, args):
    if kind == "gnn":
        e_crilin = "30.*pred_energy[0]"
    elif kind == "nognn":
        e_crilin = "tot_energy"

    e_hcal = f"VD_energy[0]*(1 + rng*{args.hcal_smear}/sqrt(VD_energy[0]))"

    return f"(({e_crilin}) + ({e_hcal}))"


def make_hist_from_tree(tree, expr, selection, name, nbins):
    tmp_name = f"tmp_{name}"
    ROOT.gDirectory.Delete(f"{tmp_name};*")
    tree.Draw(f"{expr}>>{tmp_name}(1000,-1.,1.)", selection, "goff")
    htmp = ROOT.gDirectory.Get(tmp_name)

    mean = htmp.GetMean()
    rms = htmp.GetRMS()

    ROOT.gDirectory.Delete(f"{name};*")
    h = ROOT.TH1F(name, name, int(nbins), mean-8.*rms, mean+8.*rms)
    tree.Draw(f"{expr}>>{name}", selection, "goff")
    h.SetDirectory(0)
    return h


def dcb_fit(h, name, title, root_file=None, outdir="."):

    x = ROOT.RooRealVar(f"x_{name}", "E/E_{BEAM} - 1", h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
    data = ROOT.RooDataHist(f"data_{name}", f"data_{name}", ROOT.RooArgList(x), h)

    mean = ROOT.RooRealVar(f"mean_{name}", "DCB mean", h.GetMean(), h.GetMean()-1., h.GetMean()+1.)
    sigma = ROOT.RooRealVar(f"sigma_{name}", "DCB sigma", h.GetRMS(), 0.1*h.GetRMS(), 5*h.GetRMS())
    alphaL = ROOT.RooRealVar(f"alphaL_{name}", "alphaL", 1.5, 0.1, 5.)
    nL = ROOT.RooRealVar(f"nL_{name}", "nL", 3.0, 0.5, 20.)
    alphaR = ROOT.RooRealVar(f"alphaR_{name}", "alphaR", 1.5, 0.1, 5.)
    nR = ROOT.RooRealVar(f"nR_{name}", "nR", 3.0, 0.5, 20.)

    dcb = ROOT.RooCrystalBall(f"dcb_{name}", "Double Crystal Ball", x, mean, sigma, alphaL, nL, alphaR, nR)
    nsig = ROOT.RooRealVar(f"nsig_{name}", "signal yield", h.Integral(), 0., 10.0*h.Integral())
    model = ROOT.RooAddPdf(f"model_{name}", "extended DCB model", ROOT.RooArgList(dcb), ROOT.RooArgList(nsig))

    result = model.fitTo(data, ROOT.RooFit.Extended(True), ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1))

    #FWHM
    x_min = x.getMin()
    x_max = x.getMax()
    tf1 = model.asTF(ROOT.RooArgList(x))
    halfmax = tf1.GetMaximum(x_min, x_max)/2.
    x_peak = tf1.GetMaximumX(x_min, x_max)
    x_left = tf1.GetX(halfmax, x_min, x_peak)
    x_right = tf1.GetX(halfmax, x_peak, x_max)
    fwhm_over_235 = (x_right - x_left)/2.35
    fwhm_over_235_err = abs(fwhm_over_235)*sigma.getError()/sigma.getVal()

    frame = x.frame(ROOT.RooFit.Title(title))
    data.plotOn(frame, ROOT.RooFit.MarkerSize(0.7))
    model.plotOn(frame, ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(9))
    frame.GetXaxis().SetTitle("E_{RECO}/E_{BEAM} - 1")
    frame.GetYaxis().SetTitle("Entries")
    frame.GetYaxis().SetTitleOffset(1.5)
    frame.GetXaxis().SetTitleOffset(1.2)

    c = ROOT.TCanvas(f"c_{name}", title, 800, 600)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    frame.Draw()

    pt = ROOT.TPaveText(0.58, 0.64, 0.88, 0.88, "NDC")
    pt.SetFillColor(0)
    pt.SetBorderSize(0)
    pt.SetTextFont(42)
    pt.SetTextSize(0.03)
    pt.SetTextAlign(32)
    pt.AddText(f"m_{{core}}^{{DCB}} = {mean.getVal():.3f} #pm {mean.getError():.3f}")
    pt.AddText(f"#sigma_{{core}}^{{DCB}} = {sigma.getVal():.3f} #pm {sigma.getError():.3f}")
    pt.AddText(f"FWHM / 2.35 = {fwhm_over_235:.3f} #pm {fwhm_over_235_err:.3f}")
    pt.Draw()
    c._keep = [pt, frame]
    c.Modified()
    c.Update()

    c.SaveAs(os.path.join(outdir, f"{name}.pdf"))

    if root_file:
        root_file.cd()
        h.Write(f"hist_{name}", ROOT.TObject.kOverwrite)
        c.Write(f"canvas_{name}", ROOT.TObject.kOverwrite)
        result.Write(f"fitresult_{name}", ROOT.TObject.kOverwrite)

    return {
        "mean": (mean.getVal(), mean.getError()),
        "sigma": (sigma.getVal(), sigma.getError()),
        "fwhm_over_235": (fwhm_over_235, fwhm_over_235_err),
        "fwhm_height": halfmax,
        "fwhm_xleft": x_left,
        "fwhm_xright": x_right,
        "canvas": c,
        "hist": h,
    }


def run_differential(tree, kind, args, root_file):
    x_vals, ex_vals = [], []
    bias_vals, ebias_vals = [], []
    res_vals, eres_vals = [], []

    kind_label = "With ParticleNet" if kind == "gnn" else "Without ParticleNet"
    expr = f"({reco_expr(kind, args)})/(PrimaryEnergy[0]) - 1"
    base_cut = list(args.cut)
    if args.exclude_mips:
        base_cut.append("pred_energy[0]!=0")

    print(f"\n*** Differential resolution: {kind_label} ***")

    Ebins=list(range(5, 100, 5))

    for i in range(len(Ebins) - 1):
        e_cut = f"(PrimaryEnergy[0]>{Ebins[i]} && PrimaryEnergy[0]<={Ebins[i+1]})"
        full_cut = make_cut(base_cut + [e_cut], args.verbose)
        name = f"{kind_label}_resolution_E{int(Ebins[i])}To{int(Ebins[i+1])}"
        title = f"{kind_label}: {Ebins[i]} GeV < E_{{BEAM}} < {Ebins[i+1]} GeV"

        print("Processing", name, "with cut:", full_cut)

        nbins=100
        if Ebins[i] < 20:
            nbins=25 if kind=="gnn" else 10
        if Ebins[i] < 70:
            nbins=50

        h = make_hist_from_tree(tree, expr, full_cut, f"h_{name}", nbins)
        if not h or h.GetEntries()<100:
            n = h.GetEntries() if h else 0
            print(f"** WARNING: skip {name}, entries={n}")
            continue

        result = dcb_fit(h, name, title, root_file=root_file, outdir=args.outdir)
        if result is None:
            continue

        x_vals.append((Ebins[i]+Ebins[i+1])/2.)
        ex_vals.append(0.)
        bias_vals.append(result["mean"][0])
        ebias_vals.append(result["mean"][1])
        res_vals.append(result["sigma"][0])
        eres_vals.append(result["sigma"][1])

    graphs = make_summary_plots(kind_label, x_vals, ex_vals, bias_vals, ebias_vals, res_vals, eres_vals, args, root_file)
    return graphs


def make_graph(name, title, x, ex, y, ey, marker, color):
    xa = array.array("d", x)
    exa = array.array("d", ex)
    ya = array.array("d", y)
    eya = array.array("d", ey)
    gr = ROOT.TGraphErrors(len(x), xa, ya, exa, eya)
    gr.SetName(name)
    gr.SetTitle(title)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(0.7)
    gr.SetMarkerColor(color)
    gr.SetLineColor(color)
    return gr


def make_summary_plots(kind_label, x, ex, bias, ebias, res, eres, args, root_file):
    if len(x) == 0:
        return None

    color = ROOT.kBlue+1 if kind_label == "GNN" else ROOT.kRed+1
    marker = 20

    gbias = make_graph(
        f"g_bias_{kind_label}",
        f"Bias {kind_label};E_{{BEAM}} [GeV];Bias of E_{{RECO}}/E_{{BEAM}} - 1",
        x, ex, bias, ebias, marker, color)
    gres = make_graph(
        f"g_resolution_{kind_label}",
        f"Resolution {kind_label};E_{{BEAM}} [GeV];#sigma_{{E_{{RECO}}}}/E_{{BEAM}}",
        x, ex, res, eres, marker, color)

    #bias
    cb = ROOT.TCanvas(f"c_bias_{kind_label}", f"Bias {kind_label}", 800, 600)
    cb.SetLeftMargin(0.15)
    cb.SetBottomMargin(0.15)
    gbias.Draw("AP")
    gbias.GetYaxis().SetTitleOffset(1.5)
    gbias.GetXaxis().SetTitleOffset(1.1)
    cb.Update()
    cb.SaveAs(os.path.join(args.outdir, f"bias_vs_E_{kind_label}.pdf"))

    #resol
    cr = ROOT.TCanvas(f"c_resolution_{kind_label}", f"Resolution {kind_label}", 800, 600)
    cr.SetLeftMargin(0.15)
    cr.SetBottomMargin(0.15)
    gres.Draw("AP")
    gres.GetYaxis().SetTitleOffset(1.5)
    gres.GetXaxis().SetTitleOffset(1.1)
    gres.GetYaxis().SetRangeUser(0.0, max(0.2, 1.2*max(res)))

    fit = ROOT.TF1(f"fit_res_{kind_label}", "[0]/sqrt(x)", 1., 100.)
    fit.SetParameter(0, 0.3)
    fit.SetLineWidth(2)
    fit.SetLineStyle(9)
    fit.SetLineColor(color)
    gres.Fit(fit, "R")
    fit.Draw("SAME")

    s_val = fit.GetParameter(0)
    s_err = fit.GetParError(0)
    pt = ROOT.TPaveText(0.35, 0.70, 0.88, 0.85, "NDC")
    pt.SetFillColor(0)
    pt.SetBorderSize(0)
    pt.SetTextFont(42)
    pt.SetTextSize(0.040)
    pt.SetTextAlign(32)
    color_code = 4 if kind_label == "GNN" else 2
    pt.AddText(f"#color[{color_code}]{{#font[62]{{#sigma_{{E_{{CRILIN}} + E_{{HCAL}}}} / E_{{BEAM}} = S /#sqrt{{E_{{BEAM}} [GeV]}}}}}}")
    pt.AddText(f"S [#sqrt{{GeV}}] = {s_val:.3f} #pm {s_err:.3f}")
    pt.Draw("SAME")
    cr._keep = [pt, fit]
    cr.Update()
    cr.SaveAs(os.path.join(args.outdir, f"resolution_vs_E_{kind_label}.pdf"))

    root_file.cd()
    gbias.Write(gbias.GetName(), ROOT.TObject.kOverwrite)
    gres.Write(gres.GetName(), ROOT.TObject.kOverwrite)
    cb.Write(cb.GetName(), ROOT.TObject.kOverwrite)
    cr.Write(cr.GetName(), ROOT.TObject.kOverwrite)
    fit.Write(fit.GetName(), ROOT.TObject.kOverwrite)

    return {"bias": gbias, "resolution": gres, "fit": fit, "bias_canvas": cb, "resolution_canvas": cr}


def make_overlay_en_slice(tree, args, root_file, en_center, en_width):
    base_cut = list(args.cut)
    if args.exclude_mips:
        base_cut.append(f"pred_energy[0]!=0 && abs(PrimaryEnergy[0]-{en_center})<{en_width}")
    else:
        base_cut.append(f"abs(PrimaryEnergy[0]-{en_center})<{en_width}")
    cut = make_cut(base_cut, args.verbose)

    hists = []
    for kind in ["gnn", "nognn"]:
        kind_label = "With ParticleNet" if kind == "gnn" else "Without ParticleNet"
        expr = f"({reco_expr(kind, args)})/(PrimaryEnergy[0]) - 1"
        print("\n\n\n *** Overlay 20 GeV cut:", cut)
        print("Overlay 20 GeV entries:", tree.GetEntries(cut))
        h = make_hist_from_tree(tree, expr, cut, f"h_overlay20_{kind_label}", 100)

        if not h:
            print(f"WARNING: no histogram for overlay {kind_label}")
            continue
        h.SetTitle(";(E_{CRILIN}^{RECO}+E_{VD-HCAL}^{25%/#sqrt{E}})/E_{BEAM} - 1;Entries")
        h.SetLineWidth(2)
        h.SetLineColor(ROOT.kBlue+1 if kind == "gnn" else ROOT.kRed+1)
        h.SetMarkerColor(h.GetLineColor())
        if kind=="gnn": fitres = dcb_fit(h, f"overlay20_{kind_label}", f"{kind_label}: Energy resolution (20 GeV #pi^{{+}})", root_file=root_file, outdir=args.outdir)
        else: fitres = None
        hists.append((kind_label, h, fitres))

    if len(hists) == 0:
        return None

    c = ROOT.TCanvas("c_overlay_resolution_slice", "comparison", 800, 600)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)

    maxy = max(h.GetMaximum() for _, h, _ in hists)
    leg = ROOT.TLegend(0.60, 0.70, 0.88, 0.86)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetTextFont(42)

    for i, (label, h, fitres) in enumerate(hists):

        if fitres is not None:
          #FWHM
          fwhm = fitres["fwhm_over_235"][0]
          h_to_plot = fitres["hist"]
          leg.AddEntry(h_to_plot, f"{label}: #sigma_{{DCB}} = {fitres['sigma'][0]*1e2:.1f}%, FWHM/2.35 = {1e2*fwhm:.1f}%", "l")
          l = ROOT.TLatex(fitres["fwhm_xleft"]-0.05, fitres["fwhm_height"], f"FWHM/2.35 = {fitres['sigma'][0]*1e2:.1f}%")
          l.Draw()
          a = ROOT.TArrow(fitres["fwhm_xleft"],fitres["fwhm_height"],fitres["fwhm_xright"],fitres["fwhm_height"])
          a.Draw()
        else:
          h_to_plot = h
          leg.AddEntry(h, f"{label}")

        h_to_plot.SetMaximum(1.25*maxy)
        if i==0:
            h_to_plot.SetTitle(f"Simulated events with no event cuts applied [{en_center-en_width}-{en_center+en_width} GeV #pi^{{+}}]")
            h_to_plot.Draw("HIST")
            h_to_plot.GetYaxis().SetTitleOffset(1.5)
            h_to_plot.GetXaxis().SetTitleOffset(1.1)
        else:
            h_to_plot.Draw("HIST SAME")


    leg.Draw()
    l1 = ROOT.TLatex(-0.8, 1.10*maxy, "Smearing (0.2 pe/MeV) & 40 MeV threshold on CRILIN hits")
    l1.SetTextSize(0.02)
    l1.SetTextFont(42)
    l1.Draw()
    l2  = ROOT.TLatex(-0.8, 1.05*maxy, "Virtual detector (VD) as HCAL downstream")
    l2.SetTextSize(0.02)
    l2.SetTextFont(42)
    l2.Draw()

    c.Update()
    c.SaveAs(os.path.join(args.outdir, "overlay_res_20GeV.pdf"))

    root_file.cd()
    for label, h, _ in hists:
        h.Write(f"overlay20_{label}", ROOT.TObject.kOverwrite)
    c.Write("canvas_overlay_res_20GeV", ROOT.TObject.kOverwrite)
    return c


def main():
    parser = argparse.ArgumentParser(description="Energy resolution anlysis for CRILIN GNN")
    parser.add_argument("input_file", help="input .root file")
    parser.add_argument("-o", "--output", default="gnn_reso_analysis", help="output prefix")
    parser.add_argument("--outdir", default="gnn_reso_plots", help="directory for plots")
    parser.add_argument("-c", "--cut", action="append", default=["1"], help="additional cuts")
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")

    parser.add_argument("--hcal-smear", type=float, default=0.25, help="HCAL smearing, fractional. Default: 0.25")
    parser.add_argument("--exclude-mips", action="store_true", help="apply MIP cut")

    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True) #does not open graph windows!!

    ROOT.gROOT.LoadMacro("MyRootStyle.C") #temporaneo, metti stile con una funzione
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)

    os.makedirs(args.outdir, exist_ok=True)
    root_path = f"{args.output}.root"

    tfile = ROOT.TFile.Open(args.input_file)
    tree = tfile.Get("events")

    tree.SetAlias("rng", "sin(2*pi*rndm)*sqrt(-2*log(rndm))")
    tree.SetAlias("tot_energy", "Sum$(hit_energy * (hit_energy > 0.04) * (1 + rng/sqrt(hit_energy*1e3*0.2)))")
    tree.SetAlias("hit_energy", "Hit_NCherenkov/24.3e3")

    root_file = ROOT.TFile(root_path, "RECREATE")

    #run_differential(tree, "gnn", args, root_file)
    #run_differential(tree, "nognn", args, root_file)
    make_overlay_en_slice(tree, args, root_file, 20, 5)

    root_file.Close()
    tfile.Close()
    print("*** ANALYSIS FINISHED ***")


if __name__ == "__main__":
    main()
