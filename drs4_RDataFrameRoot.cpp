// root_analysis.cpp
// Analisi del TTree "events" prodotto da drs4_ROOT.cpp con RDataFrame.
// Per ogni canale richiesto mostra (e salva in PNG) una canvas 2x2:
// ampiezza, rumore di baseline, carica e carica vs ampiezza.
//
// Uso (senza -b e senza -q, altrimenti le finestre non restano aperte):
//   root -l 'root_analysis.cpp("run.root", 0.02, "1,2")'
//   (argomenti: file root, soglia sull'ampiezza [V], canali da plottare separati da virgola)

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"

// range degli istogrammi (nelle unita' del plot: mV e pC): modificali qui
constexpr int    NB_AMP = 200;  constexpr double AMP_MIN_MV = -20.0, AMP_MAX_MV = 500.0;
constexpr int    NB_Q   = 200;  constexpr double Q_MIN_PC   = -5.0,  Q_MAX_PC   = 100.0;
constexpr int    NB_BN  = 100;  constexpr double BN_MAX_MV  = 20.0;
constexpr double RATE_BIN = 10.0;                                    // [s]

struct ChHists {
  int c;
  ROOT::RDF::RResultPtr<TH1D> amp, charge, bnoise;
  ROOT::RDF::RResultPtr<TH2D> qa;
  ROOT::RDF::RResultPtr<ROOT::RDF::RCutFlowReport> report;
};

void root_analysis(const char* fname = "run.root", double thr = 0.02,
                   const char* chans = "1,2,3,4") {
  // gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);     // le copie degli istogrammi non finiscono in gDirectory
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame df("events", fname);

  const auto cols = df.GetColumnNames();
  auto has = [&](const std::string& n) {
    return std::find(cols.begin(), cols.end(), n) != cols.end();
  };

  // tempo relativo all'inizio del run
  const double t0 = df.Min("unixtime").GetValue();
  const double t1 = df.Max("unixtime").GetValue();
  auto d = df.Define("t_rel", [t0](double t) { return t - t0; }, {"unixtime"});

  const int nbRate = std::max(1, int((t1 - t0) / RATE_BIN));
  auto hRate = d.Histo1D({"hRate", Form("Rate degli eventi;t [s];Eventi / %.0f s", RATE_BIN),
                          nbRate, 0, t1 - t0 + 1e-3}, "t_rel");
  auto nTot = d.Count();

  // canali richiesti, es. "1,2" -> {1, 2}
  std::vector<int> wanted;
  {
    std::stringstream ss(chans);
    std::string tok;
    while (std::getline(ss, tok, ','))
      if (!tok.empty()) wanted.push_back(std::stoi(tok));
  }

  // larghezze dei bin, per le etichette degli assi y
  const double wAmp = (AMP_MAX_MV - AMP_MIN_MV) / NB_AMP;
  const double wQ   = (Q_MAX_PC - Q_MIN_PC) / NB_Q;
  const double wBn  = BN_MAX_MV / NB_BN;

  // istogrammi per canale (tutti prenotati qui, un solo ciclo sugli eventi)
  std::vector<ChHists> hs;
  for (int c : wanted) {
    const std::string pre = "ch" + std::to_string(c);
    const std::string a = pre + "_amp", q = pre + "_charge", b = pre + "_bnoise";
    if (!has(a) || !has(q) || !has(b)) {
      std::cerr << "Canale " << c << " non presente nel tree: salto.\n";
      continue;
    }

    // colonne nelle unita' del plot (il tree resta in V e C)
    auto dc = d.Define(pre + "_amp_mV", a + " * 1e3")
               .Define(pre + "_q_pC",   q + " * 1e12")
               .Define(pre + "_bn_mV",  b + " * 1e3");
    const std::string aM = pre + "_amp_mV", qP = pre + "_q_pC", bM = pre + "_bn_mV";

    auto sel = dc.Filter(a + " > " + std::to_string(thr), pre + " amp > soglia");

    ChHists h;
    h.c = c;
    h.amp    = dc.Histo1D({Form("hAmp%d", c),
                           Form("Canale %d: ampiezza;Ampiezza [mV];Eventi / %.1f mV", c, wAmp),
                           NB_AMP, AMP_MIN_MV, AMP_MAX_MV}, aM);
    h.bnoise = dc.Histo1D({Form("hBn%d", c),
                           Form("Canale %d: rumore di baseline;Rumore [mV];Eventi / %.2f mV", c, wBn),
                           NB_BN, 0, BN_MAX_MV}, bM);
    h.charge = sel.Histo1D({Form("hQ%d", c),
                            Form("Canale %d: carica (amp > %.3g mV);Carica [pC];Eventi / %.1f pC",
                                 c, thr * 1e3, wQ),
                            NB_Q, Q_MIN_PC, Q_MAX_PC}, qP);
    h.qa     = sel.Histo2D({Form("hQA%d", c),
                            Form("Canale %d: carica vs ampiezza (amp > %.3g mV);Ampiezza [mV];Carica [pC]",
                                 c, thr * 1e3),
                            100, AMP_MIN_MV, AMP_MAX_MV, 100, Q_MIN_PC, Q_MAX_PC}, aM, qP);
    h.report = sel.Report();
    hs.push_back(h);
  }

  if (hs.empty()) { std::cerr << "Nessun canale trovato nel tree.\n"; return; }

  // il primo accesso ai risultati fa partire il ciclo sugli eventi
  std::cout << "Eventi totali: " << *nTot << "\n";

  // stile (da impostare prima di creare le canvas)
  gStyle->SetOptStat("emr");               // entries, media, rms
  gStyle->SetPadGridX(true);  gStyle->SetPadGridY(true);
  gStyle->SetPadTickX(1);     gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.13);   gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.13); gStyle->SetPadTopMargin(0.09);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetPalette(kViridis);

  auto axes = [](TH1& h) {
    h.GetXaxis()->SetTitleSize(0.045);  h.GetXaxis()->SetLabelSize(0.04);
    h.GetYaxis()->SetTitleSize(0.045);  h.GetYaxis()->SetLabelSize(0.04);
    h.GetXaxis()->SetTitleOffset(1.1);  h.GetYaxis()->SetTitleOffset(1.3);
  };
  auto style1d = [&](TH1& h, int fill, int line) {
    axes(h);
    h.SetFillColorAlpha(fill, 0.55);
    h.SetLineColor(line);
    h.SetLineWidth(2);
  };

  // draw: disegna nel pad CORRENTE (quello impostato da cd()) una copia dell'istogramma,
  // che sopravvive alla distruzione dei risultati di RDataFrame.
  // (DrawClone no: usa il pad "selezionato", che cd() non cambia.)
  // Le canvas sono allocate con new, altrimenti si chiuderebbero a fine funzione.
  auto draw = [](auto& res, const char* opt) {
    TH1* p = (TH1*)res->Clone();
    p->SetDirectory(nullptr);
    p->Draw(opt);
  };

  style1d(*hRate, kBrown + 1, kBrown + 3);
  TCanvas* cr = new TCanvas("rate", "rate", 900, 650);
  cr->cd();
  draw(hRate, "hist");
  cr->Update();
  cr->SaveAs("rate.png");

  for (auto& h : hs) {
    std::cout << "\n--- Canale " << h.c << " ---\n";
    h.report->Print();
    std::cout << "entries: amp " << h.amp->GetEntries()
              << ", bnoise " << h.bnoise->GetEntries()
              << ", carica " << h.charge->GetEntries() << "\n";

    style1d(*h.amp,    kAzure + 1,  kAzure + 3);
    style1d(*h.bnoise, kTeal + 1,   kTeal + 3);
    style1d(*h.charge, kOrange + 1, kOrange + 3);
    axes(*h.qa);
    h.qa->SetStats(0);                          // niente stat box sulla mappa 2D

    // una canvas 2x2 per canale: ampiezza, rumore, carica, carica vs ampiezza
    TCanvas* cv = new TCanvas(Form("canale%d", h.c), Form("Canale %d", h.c), 1400, 1000);
    cv->Divide(2, 2);

    cv->cd(1);  draw(h.amp,    "hist");
    cv->cd(2);  draw(h.bnoise, "hist");

    if (h.charge->GetEntries() > 0) {
      cv->cd(3);  draw(h.charge, "hist");
      cv->cd(4);  gPad->SetRightMargin(0.14);  gPad->SetGrid(0, 0);  draw(h.qa, "colz");
    } else {                                    // nessun evento sopra soglia
      for (int p = 3; p <= 4; ++p) {
        cv->cd(p);
        TLatex t;  t.SetNDC();  t.SetTextAlign(22);
        t.DrawLatex(0.5, 0.5, Form("nessun evento con amp > %.3g mV", thr * 1e3));
      }
    }

    cv->Update();
    cv->SaveAs(Form("ch%d_plots.png", h.c));
  }
}