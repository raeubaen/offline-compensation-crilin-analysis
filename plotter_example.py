#python3

#questo fa un dataset solo (anche molti file) ma una sola tchain

#va chiamato da un'altro programma lasciando plot.conf costante e variando data.conf per fare tutti i dataset in batch con gli stessi plot e gli stessi tagli

import ROOT
import sys
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Online monitor and reconstruction for crilin July')

parser.add_argument('plotconffile', type=str, help='plotconffile')
parser.add_argument('dataconffile', type=str, help='dataconffile')
parser.add_argument('outputfolder', type=str, help='outfolder')
parser.add_argument('--applysingleecut', type=int, help='Single particle cut', default=1)
parser.add_argument('--applycentroidcut', type=int, help='2.5mm square centroid cut', default=0)

args = parser.parse_args()
v = vars(args)
print(v)
vars().update(v)

macro = ["root_logon.C"] #x elisa: root logon

def plot(row, chain, outputfolder):
  name = row['name']
  print(f"{outputfolder}/{name}.root")
  f = ROOT.TFile(f"{outputfolder}/{name}.root", "recreate")
  f.cd()
  c = ROOT.TCanvas(f"{name}_canvas")
  c.cd()
  if str(row.cuts) == "": cut = "1"
  else: cut = str(row.cuts)
  if applysingleecut: cut = cut + " && single_e_flag[0]==1"
  if applycentroidcut: cut = cut + " && centroid_cut_flag[0]==1"
  if str(row.y).strip()=="0":
    h = ROOT.TH1F(f"{name}", f"{row.title}", int(row.binsnx), float(row.binsminx), float(row.binsmaxx))
    chain.Draw(f"{row.x}>>{name}", f"{cut}")
    print(f"{row.x}>>{name}, {cut}")
    h.SetLineColor(eval(f"ROOT.{row.color}"))
    binw = (float(row.binsmaxx) - float(row.binsminx))/int(row.binsnx)
    h.GetYaxis().SetTitle(f"Entries / {float(f'{binw:.1g}'):g} {row.ylabel}")
  else:
    h = ROOT.TH2F(f"{name}", f"{row.title}", int(row.binsnx), float(row.binsminx), float(row.binsmaxx), int(row.binsny), float(row.binsminy), float(row.binsmaxy))
    chain.Draw(f"{row.y}:{row.x}>>{name}", f"{cut}", "zcol")
    h.GetYaxis().SetTitle(f"{row.ylabel}")
  h.GetXaxis().SetTitle(f"{row.xlabel}")
  c.SaveAs(f"{outputfolder}/{name}.pdf")
  c.SaveAs(f"{outputfolder}/{name}.png")
  c.Write()
  h.Write()
  f.Close()
  c.Close()
  del c
  del h

def process(row, outputfolder, plot_df):
  lst = os.popen(f"/bin/bash -c 'ls -1 {row.filename.strip()}'").read().split("\n")
  chain = ROOT.TChain()
  for file in lst:
    chain.Add(f"{file}/{row.treename.strip()}")
  os.system(f"mkdir {outputfolder}/{row.label}")
  os.system(f"cp index.php {outputfolder}/{row.label}")
  plotconf_df.apply(lambda plotrow: plot(plotrow, chain, f"{outputfolder}/{row.label}"), axis=1)

dataconf_df = pd.read_csv(dataconffile, sep=";")

for m in macro: ROOT.gROOT.LoadMacro(m)

plotconf_df = pd.read_csv(plotconffile, sep=";")

os.system(f"cp index.php {outputfolder}")

dataconf_df.apply(lambda row: process(row, outputfolder, plotconf_df), axis=1)
