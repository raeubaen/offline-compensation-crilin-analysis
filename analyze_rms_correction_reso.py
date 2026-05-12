import argparse
import array
import os
import ROOT


def make_cut(cuts):
    cuts = ["("+c+")" for c in cuts]
    cut="*".join(cuts)
    return cut


def setup_aliases(tree, hit_energy_cut):
    cherenkov_scale = 24.3e3

    tree.SetAlias("rng", "sin(2*pi*rndm)*sqrt(-2*log(rndm))")
    tree.SetAlias("hit_energy", f"Hit_NCherenkov/{cherenkov_scale}")

    tree.SetAlias("variance", f"Sum$((Hit_x*Hit_x + Hit_y*Hit_y)*Hit_NCherenkov*(hit_energy>{hit_energy_cut}))/Sum$(Hit_NCherenkov * (hit_energy>{hit_energy_cut}))")
    tree.SetAlias("RMS_cluster", "sqrt(variance)")

    tree.SetAlias("CRILIN_true", "PrimaryEnergy[0] - VD_energy[0]")
    tree.SetAlias("CRILIN_cherenkov_calib", f"Sum$(Hit_NCherenkov)/{cherenkov_scale}")
    tree.SetAlias("CRILIN_response", "CRILIN_cherenkov_calib / CRILIN_true")

    tree.SetAlias("CRILIN_cherenkov_photostat", f"Sum$(hit_energy * (hit_energy > {hit_energy_cut}) * (1 + rng/sqrt(hit_energy*1e3*0.2)))")


def make_hist_from_tree(tree, expr, selection, name, nbins, scan_range=(-1.0, 1.0), nsigma_range=6.0):
    tmp_name = f"tmp_{name}"
    ROOT.gDirectory.Delete(f"{tmp_name};*")
    ROOT.gDirectory.Delete(f"{name};*")

    n_tmp = tree.Draw(f"{expr}>>{tmp_name}(1000,{scan_range[0]},{scan_range[1]})", selection, "goff")
    htmp = ROOT.gDirectory.Get(tmp_name)

    if n_tmp <= 0 or not htmp or not htmp.InheritsFrom("TH1") or htmp.GetEntries() < 10:
        print(f"WARNING: cannot build temporary histogram {tmp_name}")
        print("  expr =", expr)
        print("  cut  =", selection)
        return None

    mean = htmp.GetMean()
    rms = htmp.GetRMS()
    if rms <= 0:
        print(f"WARNING: RMS <= 0 for {tmp_name}")
        return None

    h = ROOT.TH1F(name, ";(E_{CRILIN}^{corr} + E_{HCAL})/E_{BEAM} - 1;Entries", int(nbins), mean - nsigma_range * rms, mean + nsigma_range * rms)
    h.Sumw2()

    n = tree.Draw(f"{expr}>>{name}", selection, "goff")
    if n <= 0 or h.GetEntries() <= 0:
        print(f"WARNING: final histogram {name} is empty")
        return None

    h.SetDirectory(0)
    return h


def dcb_fit(h, name, title, root_file=None, outdir="."):
    if not h or h.GetEntries() < 20:
        print(f"WARNING: dcb_fit skipped for {name}")
        return None

    x = ROOT.RooRealVar(f"x_{name}",  "E/E_{BEAM} - 1", h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
    data = ROOT.RooDataHist(f"data_{name}", f"data_{name}", ROOT.RooArgList(x), h)

    h_mean = h.GetMean()
    h_rms = max(h.GetRMS(), 1e-6)

    mean = ROOT.RooRealVar(f"mean_{name}", "DCB mean", h_mean, h_mean - 1.0, h_mean + 1.0)
    sigma = ROOT.RooRealVar(f"sigma_{name}", "DCB sigma", h_rms, 0.05 * h_rms, 10.0 * h_rms)
    alphaL = ROOT.RooRealVar(f"alphaL_{name}", "alphaL", 1.5, 0.1, 8.0)
    nL = ROOT.RooRealVar(f"nL_{name}", "nL", 3.0, 0.5, 50.0)
    alphaR = ROOT.RooRealVar(f"alphaR_{name}", "alphaR", 1.5, 0.1, 8.0)
    nR = ROOT.RooRealVar(f"nR_{name}", "nR", 3.0, 0.5, 50.0)

    dcb = ROOT.RooCrystalBall(f"dcb_{name}", "Double Crystal Ball", x, mean, sigma, alphaL, nL, alphaR, nR)
    nsig = ROOT.RooRealVar(f"nsig_{name}", "signal yield", h.Integral(), 0., 10.*h.Integral())
    model = ROOT.RooAddPdf(f"model_{name}", "extended DCB model", ROOT.RooArgList(dcb), ROOT.RooArgList(nsig))

    obs = ROOT.RooArgSet(x)
    model.fixCoefNormalization(obs)

    result = model.fitTo(data, ROOT.RooFit.Extended(True), ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1))

    #FWHM

    probs = array.array("d", [0.16, 0.84])
    quants = array.array("d", [0,0])

    h.GetQuantiles(2, quants, probs)

    fwhm_over_235 = (quants[1]-quants[0])/2.

    if sigma.getVal()>0:
        fwhm_over_235_err = abs(fwhm_over_235)*sigma.getError()/sigma.getVal()
    else:
        fwhm_over_235_err = 0.

    frame = x.frame(ROOT.RooFit.Title(title))
    data.plotOn(frame, ROOT.RooFit.MarkerSize(0.7))
    model.plotOn(frame, ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(9))
    frame.GetXaxis().SetTitle("(E_{CRILIN}^{corr}+E_{HCAL})/E_{BEAM} - 1")
    frame.GetYaxis().SetTitle("Entries")

    c = ROOT.TCanvas(f"c_{name}", title, 800, 600)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    frame.Draw()

    pt = ROOT.TPaveText(0.55, 0.64, 0.88, 0.88, "NDC")
    pt.SetFillColor(0)
    pt.SetBorderSize(0)
    pt.SetTextFont(42)
    pt.SetTextSize(0.03)
    pt.SetTextAlign(32)
    pt.AddText(f"m_{{core}}^{{DCB}} = {mean.getVal():.4f} #pm {mean.getError():.4f}")
    pt.AddText(f"#sigma_{{core}}^{{DCB}} = {sigma.getVal():.4f} #pm {sigma.getError():.4f}")
    pt.AddText(f"FWHM / 2.35 = {fwhm_over_235:.4f} #pm {fwhm_over_235_err:.4f}")
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

    return_dict = {
        "mean": (mean.getVal(), mean.getError()),
        "sigma": (sigma.getVal(), sigma.getError()),
        "fwhm_over_235": (fwhm_over_235, fwhm_over_235_err),
        "fit_status": result.status(),
        "cov_qual": result.covQual(),
        "canvas": c,
        "hist": h,
        "frame": frame,
        "model": model,
    }

    return return_dict

def make_resolution_comparison_plot(
        h_uncorr,
        h_corr,
        fit_corr,
        name,
        title,
        outdir,
        root_file=None):

    c = ROOT.TCanvas(f"c_compare_{name}", title, 850, 700)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.13)

    maxy = max(h_uncorr.GetMaximum(), h_corr.GetMaximum()) * 1.35

    h_uncorr.SetLineColor(ROOT.kBlack)
    h_uncorr.SetLineWidth(2)
    h_uncorr.SetMarkerColor(ROOT.kBlack)

    h_corr.SetLineColor(ROOT.kBlue + 1)
    h_corr.SetLineWidth(2)
    h_corr.SetMarkerColor(ROOT.kBlue + 1)

    h_uncorr.SetMaximum(maxy)

    h_uncorr.Draw("HIST")
    h_corr.Draw("HIST SAME")

    leg = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    leg.AddEntry(h_uncorr, "Without RMS correction", "l")
    leg.AddEntry(h_corr, "With RMS correction", "l")

    if fit_corr is not None:
        leg.AddEntry(0, f"#sigma = {fit_corr['fwhm_over_235'][0]:.4f}", "")

    leg.Draw()

    c._keep = [leg]

    c.Modified()
    c.Update()

    c.SaveAs(os.path.join(outdir, f"{name}.pdf"))

    if root_file:
        root_file.cd()
        c.Write(c.GetName(), ROOT.TObject.kOverwrite)

    return c

def make_graph(name, title, x, ex, y, ey, marker=20, color=ROOT.kBlue + 1):
    xa = array.array("d", x)
    exa = array.array("d", ex)
    ya = array.array("d", y)
    eya = array.array("d", ey)
    gr = ROOT.TGraphErrors(len(x), xa, ya, exa, eya)
    gr.SetName(name)
    gr.SetTitle(title)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(0.8)
    gr.SetMarkerColor(color)
    gr.SetLineColor(color)
    return gr


def fit_rms_response_profile(tree, cut, tag, title, args, root_file):
    h2_name = f"h2_rms_vs_response_{tag}"
    ROOT.gDirectory.Delete(f"{h2_name};*")

    xbins = 4*args.nbins_x if tag=="overall" else args.nbins_x
    ybins = 4*args.nbins_y if tag=="overall" else args.nbins_y

    h2 = ROOT.TH2F(h2_name, ";CRILIN cluster RMS [mm];E_{CRILIN}^{reco}/E_{CRILIN}^{true}", xbins, args.xmin, args.xmax, ybins, args.ymin, args.ymax)

    tree.Draw(f"CRILIN_response:RMS_cluster>>{h2_name}", cut, "goff")
    if h2.GetEntries() < args.min_entries:
        print(f"WARNING: skip RMS response fit for {tag}, entries={h2.GetEntries()}")
        return None

    prof = h2.ProfileX(f"prof_rms_vs_response_{tag}")
    prof.SetMarkerStyle(20)
    prof.SetMarkerSize(0.8)
    prof.SetMarkerColor(ROOT.kBlack)
    prof.SetLineColor(ROOT.kBlack)
    prof.SetLineWidth(2)

    fit = ROOT.TF1(f"fit_response_vs_rms_{tag}", "pol1", args.fit_xmin, args.fit_xmax)
    fit.SetLineWidth(3)
    fit.SetLineStyle(9)
    fit.SetLineColor(ROOT.kRed)

    prof.Fit(fit, "R")

    p0 = fit.GetParameter(0)
    p0_err = fit.GetParError(0)
    p1 = fit.GetParameter(1)
    p1_err = fit.GetParError(1)

    mean_rms = h2.GetMean(1)
    mean_rms_err = h2.GetMeanError(1)

    print(title)
    c = ROOT.TCanvas(f"c_rms_vs_response_{tag}", title, 900, 700)
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.15)
    c.SetBottomMargin(0.13)

    h2.Draw("COLZ")
    prof.Draw("SAME")
    fit.Draw("SAME")

    pt = ROOT.TPaveText(0.43, 0.68, 0.84, 0.88, "NDC")
    pt.SetFillColor(0)
    pt.SetBorderSize(0)
    pt.SetTextFont(42)
    pt.SetTextSize(0.033)
    pt.SetTextAlign(32)
    pt.AddText(f"#font[62]{{{title}}}")
    pt.AddText("#color[2]]{Fit response = p_{0} + p_{1} * RMS}")
    pt.AddText(f"p0 = {p0:.5g} #pm {p0_err:.2g}")
    pt.AddText(f"p1 = {p1:.4g} #pm {p1_err:.2g}")
    pt.Draw("SAME")

    c._keep = [pt, prof, fit]
    c.Modified()
    c.Update()
    c.SaveAs(os.path.join(args.outdir, f"rms_response_{tag}.pdf"))

    root_file.cd()
    h2.Write(h2.GetName(), ROOT.TObject.kOverwrite)
    prof.Write(prof.GetName(), ROOT.TObject.kOverwrite)
    fit.Write(fit.GetName(), ROOT.TObject.kOverwrite)
    c.Write(c.GetName(), ROOT.TObject.kOverwrite)

    return {
        "p0": (p0, p0_err),
        "p1": (p1, p1_err),
        "mean_rms": (mean_rms, mean_rms_err),
        "h2": h2,
        "profile": prof,
        "fit": fit,
        "canvas": c,
    }


def make_summary_plot(name, graph, args, root_file, fit_stochastic=False):
    c = ROOT.TCanvas(f"c_{name}", name, 800, 600)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.13)

    graph.Draw("AP")

    fit = None
    pt = None
    if fit_stochastic:
        yvals = [graph.GetPointY(i) for i in range(graph.GetN())]
        graph.GetYaxis().SetRangeUser(0.0, max(0.2, 1.3 * max(yvals)))

        fit = ROOT.TF1(f"fit_{name}", "sqrt([0]*[0]/x + [1]*[1]/(x*x))", max(args.energy_min, 1.0), args.energy_max)
        fit.SetParName(0, "S")
        fit.SetParameter(0, 0.836486)
        fit.SetLineWidth(2)
        fit.SetLineStyle(9)
        fit.SetLineColor(ROOT.kBlue)
        fit.SetParLimits(0, 0, 1000);
        fit.SetParLimits(1, 0, 1000);

        graph.Fit(fit, "REMQ")
        fit.Draw("SAME")

        S = fit.GetParameter(0)
        eS = fit.GetParError(0)
        pt = ROOT.TPaveText(0.45, 0.65, 0.88, 0.86, "NDC")
        pt.SetFillColor(0)
        pt.SetBorderSize(0)
        pt.SetTextFont(42)
        pt.SetTextSize(0.03)
        pt.SetTextAlign(32)
        pt.AddText("#font[62]{Resolution using RMS-based correction}")
        pt.AddText("#font[52]{MIPs excluded, 25%/#sqrt{E} HCAL smearing, 0.2 p.e./MeV}")
        pt.AddText("#font[62]{#color[4]{#sigma_{E_{reco}}/E_{BEAM} = S/#sqrt{E_{BEAM} [GeV]}}}")
        pt.AddText(f"S [#sqrt{{GeV}}] = {S:.3f} #pm {eS:.3f}")
        pt.Draw("SAME")
        c._keep = [fit, pt]

    c.Modified()
    c.Update()
    c.SaveAs(os.path.join(args.outdir, f"{name}.pdf"))

    root_file.cd()
    graph.Write(graph.GetName(), ROOT.TObject.kOverwrite)
    c.Write(c.GetName(), ROOT.TObject.kOverwrite)
    if fit:
        fit.Write(fit.GetName(), ROOT.TObject.kOverwrite)

    return c

def make_resolution_summary_comparison(
        g_uncorr,
        g_corr,
        args,
        root_file):

    c = ROOT.TCanvas("c_resolution_comparison", "", 850, 700)

    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.13)

    g_uncorr.SetMarkerStyle(24)
    g_uncorr.SetMarkerColor(ROOT.kRed)
    g_uncorr.SetLineColor(ROOT.kRed)

    g_corr.SetMarkerStyle(20)
    g_corr.SetMarkerColor(ROOT.kBlue+1)
    g_corr.SetLineColor(ROOT.kBlue+1)

    ymax = 0.
    for i in range(g_uncorr.GetN()):
        ymax = max(ymax, g_uncorr.GetPointY(i))

    for i in range(g_corr.GetN()):
        ymax = max(ymax, g_corr.GetPointY(i))

    g_uncorr.GetYaxis().SetRangeUser(0., ymax * 1.4)

    g_uncorr.Draw("AP")
    g_corr.Draw("P SAME")


    pt = ROOT.TPaveText(0.45, 0.65, 0.88, 0.86, "NDC")
    pt.SetFillColor(0)
    pt.SetBorderSize(0)
    pt.SetTextFont(42)
    pt.SetTextSize(0.03)
    pt.SetTextAlign(32)
    pt.AddText("#font[62]{Resolution with/without RMS-based correction}")
    pt.AddText("#font[52]{No cuts, 25%/#sqrt{E} HCAL smearing, 0.2 p.e./MeV}")
    pt.AddText("#font[62]{#color[4]{#sigma_{E_{reco}}/E_{BEAM} = S/#sqrt{E_{BEAM}}} #oplus N/E_{BEAM}}")

    d = {"With correction": g_corr, "Without correction": g_uncorr}
    for key in d:
      fit = ROOT.TF1(
          f"fit_resolution_{key}",
          "sqrt([0]*[0]/x + [1]*[1]/(x*x))",
          max(args.energy_min, 1.0),
          args.energy_max
      )

      fit.SetParameter(0, 0.666)
      fit.SetLineStyle(9)
      fit.SetLineWidth(3)

      fit.SetParLimits(0, 0, 1000);
      fit.SetParLimits(1, 0, 1000);

      d[key].Fit(fit, "REMQ")
      fit.Draw("SAME")

      S = fit.GetParameter(0)
      eS = fit.GetParError(0)

      print("stochasticx term: ", S)

      pt.AddText(f"{key}: S [#sqrt{{GeV}}] = {S*1e2:.1f}%")

    pt.Draw("SAME")

    leg = ROOT.TLegend(0.45, 0.65, 0.55, 0.70)

    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    leg.AddEntry(g_uncorr, "Without RMS correction", "pl")
    leg.AddEntry(g_corr, "With RMS correction", "pl")

    leg.Draw()

    c._keep = [leg, fit]

    c.Modified()
    c.Update()

    c.SaveAs(os.path.join(args.outdir, "resolution_vs_energy_comparison.pdf"))

    root_file.cd()

    c.Write(c.GetName(), ROOT.TObject.kOverwrite)
    fit.Write(fit.GetName(), ROOT.TObject.kOverwrite)

def run_analysis(tree, args, root_file, hit_energy_cut):
    base_cut = list(args.cut)
    #base_cut.append("pred_energy[0]!=0") #always exclude mips here

    #base_cut.append("CRILIN_true > 0")
    #base_cut.append(f"Sum$(Hit_NCherenkov * (hit_energy > {hit_energy_cut})) > 500")

    hcal_expr = f"VD_energy[0]*(1 + rng*{args.hcal_smear}/sqrt(VD_energy[0]))"
    crilin_uncorr_for_resolution = "CRILIN_cherenkov_photostat" if args.use_crilin_photostat else "CRILIN_cherenkov_calib"

    e_bins = []
    e = args.energy_min
    while e < args.energy_max:
        e_bins.append((e, min(e + args.energy_step, args.energy_max)))
        e += args.energy_step

    uncorr_res_vals = []
    euncorr_res_vals = []

    x_vals, ex_vals = [], []
    p0_vals, ep0_vals = [], []
    p1_vals, ep1_vals = [], []
    mean_rms_vals, emean_rms_vals = [], []
    bias_vals, ebias_vals = [], []
    res_vals, eres_vals = [], []

    #overall RMS-response plot, using the full energy range
    overall_cut = make_cut(base_cut + [f"PrimaryEnergy[0]>{args.energy_min} && PrimaryEnergy[0]<={args.energy_max}"])
    fit_rms_response_profile(tree, overall_cut, "overall", "", args, root_file)

    for emin, emax in e_bins:
        tag = f"E{int(emin)}To{int(emax)}"
        e_cut = f"PrimaryEnergy[0]>{emin} && PrimaryEnergy[0]<={emax}"
        full_cut = make_cut(base_cut + [e_cut])

        print(f"Processing {tag}")

        if tree.GetEntries(full_cut) < args.min_entries:
            print(f"WARNING: skip {tag}, too few entries")
            continue

        fit_info = fit_rms_response_profile(tree, full_cut, tag, f"{emin:.0f} < E_{{BEAM}} #leq {emax:.0f} GeV", args, root_file)
        if fit_info is None:
            continue

        p0, p0_err = fit_info["p0"]
        p1, p1_err = fit_info["p1"]
        mean_rms, mean_rms_err = fit_info["mean_rms"]

        correction_den = f"({p0:.12g} + ({p1:.12g})*RMS_cluster)"
        crilin_corr = f"(({crilin_uncorr_for_resolution})/{correction_den})"

        response_uncorr = f"(({crilin_uncorr_for_resolution}) + ({hcal_expr})) / PrimaryEnergy[0] - 1"

        response_corr = f"(({crilin_corr}) + ({hcal_expr})) / PrimaryEnergy[0] - 1"

        nbins = args.res_nbins
        if emin < 20:
            nbins = max(20, int(args.res_nbins / 2))
        if emin < 10:
            nbins = max(15, int(args.res_nbins / 3))

        h_uncorr = make_hist_from_tree(
            tree,
            response_uncorr,
            full_cut,
            f"h_uncorrected_resolution_{tag}",
            nbins,
            scan_range=(-1.0, 1.0),
            nsigma_range=6.
        )

        if not h_uncorr or h_uncorr.GetEntries() < args.min_entries:
            print(f"WARNING: skip uncorrected histogram for {tag}")
            continue

        probs = array.array("d", [0.16, 0.84])
        quants = array.array("d", [0,0])

        h_uncorr.GetQuantiles(2, quants, probs)

        uncorr_res = (quants[1] - quants[0]) / 2.

        uncorr_res_vals.append(uncorr_res)
        euncorr_res_vals.append(uncorr_res*0.05)

        h_res = make_hist_from_tree(tree, response_corr, full_cut, f"h_corrected_resolution_{tag}", nbins, scan_range=(-1.0, 1.0), nsigma_range=6.)
        if not h_res or h_res.GetEntries() < args.min_entries:
            n = h_res.GetEntries() if h_res else 0
            print(f"WARNING: skip corrected resolution fit for {tag}, entries={n}")
            continue

        res_fit = dcb_fit(h_res, f"corrected_resolution_{tag}", f"Corrected resolution: {emin:.0f} < E_{{BEAM}} #leq {emax:.0f} GeV", root_file=root_file, outdir=args.outdir)
        if res_fit is None:
            continue

        fwhm235, efwhm235 = res_fit["fwhm_over_235"]
        bias, ebias = res_fit["mean"]

        print(emin, emax,  fwhm235)
        res_vals.append(fwhm235)
        eres_vals.append(efwhm235)


        make_resolution_comparison_plot(
            h_uncorr,
            h_res,
            res_fit,
            f"resolution_comparison_{tag}",
            f"{emin:.0f} < E_{{BEAM}} #leq {emax:.0f} GeV",
            args.outdir,
            root_file
        )


        x_vals.append(0.5 * (emin + emax))
        ex_vals.append(0.)
        p0_vals.append(p0)
        ep0_vals.append(p0_err)
        p1_vals.append(p1)
        ep1_vals.append(p1_err)
        mean_rms_vals.append(mean_rms)
        emean_rms_vals.append(mean_rms_err)
        bias_vals.append(bias)
        ebias_vals.append(ebias)

        print(f"Corrected resolution {tag}: FWHM/2.35 = {fwhm235:.6g} +/- {efwhm235:.3g}")

    if len(x_vals) == 0:
        print("WARNING: no valid slices found")
        return

    g_p0_e = make_graph("g_p0_vs_energy", ";E_{BEAM} [GeV];p0", x_vals, ex_vals, p0_vals, ep0_vals, color=ROOT.kBlue)
    g_p1_e = make_graph("g_p1_vs_energy", ";E_{BEAM} [GeV];p1", x_vals, ex_vals, p1_vals, ep1_vals, color=ROOT.kRed)
    g_p0_rms = make_graph("g_p0_vs_mean_rms", ";Mean RMS_{cluster} [mm];p0", mean_rms_vals, emean_rms_vals, p0_vals, ep0_vals, color=ROOT.kBlue)
    g_p1_rms = make_graph("g_p1_vs_mean_rms", ";Mean RMS_{cluster} [mm];p1", mean_rms_vals, emean_rms_vals, p1_vals, ep1_vals, color=ROOT.kRed)
    g_bias = make_graph("g_corrected_bias_vs_energy",  ";E_{BEAM} [GeV];Bias of corrected response", x_vals, ex_vals, bias_vals, ebias_vals, color=ROOT.kBlack)

    g_res_corr = make_graph(
        "g_corrected_resolution_vs_energy",
        ";E_{BEAM} [GeV];#sigma_{E_{reco}} / E_{BEAM}",
        x_vals,
        ex_vals,
        res_vals,
        eres_vals,
        color=ROOT.kBlue+1
    )

    g_res_uncorr = make_graph(
        "g_uncorrected_resolution_vs_energy",
        ";E_{BEAM} [GeV];#sigma_{E_{reco}} / E_{BEAM}",
        x_vals,
        ex_vals,
        uncorr_res_vals,
        euncorr_res_vals,
        color=ROOT.kBlack
    )

    make_summary_plot("p0_vs_energy", g_p0_e, args, root_file)
    make_summary_plot("p1_vs_energy", g_p1_e, args, root_file)
    make_summary_plot("p0_vs_mean_rms", g_p0_rms, args, root_file)
    make_summary_plot("p1_vs_mean_rms", g_p1_rms, args, root_file)
    make_summary_plot("corrected_bias_vs_energy", g_bias, args, root_file)

    make_resolution_summary_comparison(
        g_res_uncorr,
        g_res_corr,
        args,
        root_file
    )


def main():
    parser = argparse.ArgumentParser(description="RMS-based CRILIN correction and energy resolution")

    parser.add_argument("input_file", help="input .root file containing tree 'events'")
    parser.add_argument("-o", "--output", default="rms_corrected_reso", help="output ROOT prefix")
    parser.add_argument("--outdir", default="rms_corrected_reso_plots", help="directory for PDF plots")
    parser.add_argument("-c", "--cut", action="append", default=["1"], help="additional cuts")

    parser.add_argument("--hcal-smear", type=float, default=0.25, help="HCAL smearing term")
    parser.add_argument("--use-crilin-photostat", action="store_true", default=True, help="use CRILIN photostatistical smearing in final resolution")

    parser.add_argument("--energy-min", type=float, default=5.)
    parser.add_argument("--energy-max", type=float, default=100.)
    parser.add_argument("--energy-step", type=float, default=5.)

    parser.add_argument("--nbins-x", type=int, default=40)
    parser.add_argument("--xmin", type=float, default=1.)
    parser.add_argument("--xmax", type=float, default=30.)
    parser.add_argument("--nbins-y", type=int, default=40)
    parser.add_argument("--ymin", type=float, default=0.)
    parser.add_argument("--ymax", type=float, default=1.5)
    parser.add_argument("--fit-xmin", type=float, default=7.)
    parser.add_argument("--fit-xmax", type=float, default=22.)

    parser.add_argument("--res-nbins", type=int, default=80)
    parser.add_argument("--min-entries", type=int, default=100)
    parser.add_argument("--batch", action="store_true", help="run ROOT in batch mode")

    args = parser.parse_args()

    if args.batch:
        ROOT.gROOT.SetBatch(True)

    ###############
    # MyRootStyle #
    ###############

    ROOT.gStyle.SetLabelSize(0.04, "X")
    ROOT.gStyle.SetLabelSize(0.04, "Y")
    ROOT.gStyle.SetLabelSize(0.04, "Z")
    ROOT.gStyle.SetTitleSize(0.045, "X")
    ROOT.gStyle.SetTitleSize(0.045, "Y")
    ROOT.gStyle.SetTitleSize(0.045, "Y")
    ROOT.gStyle.SetTitleFont(62, "XYZ")
    ROOT.gStyle.SetTitleFont(62, "")
    ROOT.gStyle.SetPadBottomMargin(0.15)
    ROOT.gStyle.SetPadTopMargin(0.1)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.12)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetTitleOffset(1.2, "X")
    ROOT.gStyle.SetTitleOffset(1.5, "Y")
    ROOT.gStyle.SetPalette(ROOT.kLightTemperature)
    #TColor::InvertPalette()
    ROOT.gStyle.SetNumberContours(255)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)

    os.makedirs(args.outdir, exist_ok=True)
    root_path = f"{args.output}.root"

    input_file = ROOT.TFile.Open(args.input_file)
    if not input_file or input_file.IsZombie():
        raise RuntimeError(f"Cannot open input file: {args.input_file}")

    tree = input_file.Get("events")
    if not tree:
        raise RuntimeError("Tree 'events' not found")

    hit_energy_cut = 0.04
    setup_aliases(tree, hit_energy_cut)

    root_file = ROOT.TFile(root_path, "RECREATE")
    run_analysis(tree, args, root_file, hit_energy_cut)
    root_file.Close()
    input_file.Close()

    print("*** ANALYSIS FINISHED ***")
    print("ROOT output:", root_path)
    print("PDF output directory:", args.outdir)


if __name__ == "__main__":
    main()

