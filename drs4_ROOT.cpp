// drs4_ROOT.cpp
// Legge un file binario DRS4 (formato "DRS2"), scrive un TTree con una entry
// per evento, poi mostra (e salva in PNG) gli istogrammi di ampiezza e carica
// di un canale, ciascuno in una canvas separata.
//
// Uso:
//   root -l 'drs4_ROOT.cpp("run.dat", 1, "run.root")'
//   (argomenti: file binario, canale da plottare, file root di output)

#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "TROOT.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"

constexpr int NCELL = 1024;
// finestre in ns, sull'asse dei tempi della waveform (quello dei plot)
constexpr double AMP_LO_NS = 56.0, AMP_HI_NS = 70.0;   // ricerca del minimo (ampiezza)
constexpr double BN_LO_NS  = 2.0,  BN_HI_NS  = 54.0;   // baseline
constexpr double Q_PRE_NS  = 6.0,  Q_POST_NS = 10.0;   // carica: ns prima/dopo il minimo

// istogrammi di output (nelle unita' del plot: mV e pC)
constexpr int    NB_AMP = 200;  constexpr double AMP_MIN_MV = -20.0, AMP_MAX_MV = 500.0;
constexpr int    NB_Q   = 200;  constexpr double Q_MIN_PC   = -5.0,  Q_MAX_PC   = 100.0;
constexpr bool   LOGY   = false;                       // true: asse y logaritmico

struct DrsEvent {
  uint32_t id = 0;
  uint16_t date[7] = {0};                    // y, m, d, h, min, s, ms
  uint16_t range_center = 0;
  uint16_t trig_cell = 0;
  std::map<int, uint32_t> scaler;
  std::map<int, std::vector<uint16_t>> adc;  // canale -> 1024 campioni
};

class DrsReader {
public:
  std::vector<int> channels;
  std::map<int, std::vector<float>> tw;      // time widths [ns]
  uint16_t board_id = 0;

  bool open(const char* fname) {
    f.open(fname, std::ios::binary);
    if (!f) { std::cerr << "Impossibile aprire " << fname << "\n"; return false; }
    char h[4];
    f.read(h, 4);
    if (std::strncmp(h, "DRS2", 4)) { std::cerr << "Non e' un file DRS4\n"; return false; }
    f.read(h, 4);
    if (std::strncmp(h, "TIME", 4)) { std::cerr << "Manca header TIME\n"; return false; }

    f.read(h, 4);
    if (h[0] == 'B' && h[1] == '#') {
      std::memcpy(&board_id, h + 2, 2);
      f.read(h, 4);
      while (h[0] == 'C') {
        int ch = std::atoi(std::string(h + 1, 3).c_str());
        channels.push_back(ch);
        tw[ch].resize(NCELL);
        f.read(reinterpret_cast<char*>(tw[ch].data()), NCELL * 4);
        f.read(h, 4);
      }
    }
    f.seekg(-4, std::ios::cur);              // torna all'inizio di "EHDR"
    return bool(f);
  }

  bool next(DrsEvent& ev) {
    char h[4];
    if (!f.read(h, 4)) return false;
    if (std::strncmp(h, "EHDR", 4)) return false;
    f.read(reinterpret_cast<char*>(&ev.id), 4);
    f.read(reinterpret_cast<char*>(ev.date), 14);
    f.read(reinterpret_cast<char*>(&ev.range_center), 2);
    f.read(h, 2);                                            // "B#"
    uint16_t bid; f.read(reinterpret_cast<char*>(&bid), 2);
    f.read(h, 2);                                            // "T#"
    f.read(reinterpret_cast<char*>(&ev.trig_cell), 2);
    for (int ch : channels) {
      f.read(h, 4);                                          // "Cxxx"
      f.read(reinterpret_cast<char*>(&ev.scaler[ch]), 4);
      ev.adc[ch].resize(NCELL);
      f.read(reinterpret_cast<char*>(ev.adc[ch].data()), NCELL * 2);
    }
    return bool(f);
  }

  // asse dei tempi [ns]: cumsum delle time width, ruotate sulla trigger cell
  std::vector<double> time(int ch, uint16_t tcell) const {
    std::vector<double> t(NCELL);
    double s = 0;
    for (int j = 0; j < NCELL; ++j) {
      s += tw.at(ch)[(j - 1 + tcell + NCELL) % NCELL];
      t[j] = s;
    }
    return t;
  }
private:
  std::ifstream f;
};

void drs4_ROOT(const char* fname = "run.dat", int ch = 1,
               const char* outName = "run.root") {
  // gROOT->SetBatch(kTRUE);

  DrsReader r;
  if (!r.open(fname)) return;

  const size_t nch = r.channels.size();
  std::cout << "Scheda " << r.board_id << ", canali:";
  for (int c : r.channels) std::cout << " " << c;
  std::cout << "\n";

  if (r.channels.empty()) { std::cerr << "Nessun canale nel file.\n"; return; }

  // passo di campionamento medio [ns], dalle time width di tutti i canali
  double sumTw = 0;
  for (int c : r.channels) for (float w : r.tw[c]) sumTw += w;
  // const double dt = sumTw / (double(NCELL) * r.channels.size());
  const double dt = 0.2;  // passo nominale DRS4, per convertire le finestre in ns -> campioni

  // ns -> indice di campione (limitato a [0, NCELL])
  auto ns2s = [&](double ns) { return std::min(NCELL, std::max(0, (int)std::lround(ns / dt))); };
  const int AMP_LO = ns2s(AMP_LO_NS), AMP_HI = ns2s(AMP_HI_NS);
  const int BN_LO  = ns2s(BN_LO_NS),  BN_HI  = ns2s(BN_HI_NS);
  const int Q_PRE  = ns2s(Q_PRE_NS),  Q_POST = ns2s(Q_POST_NS);   // larghezze in campioni

  std::cout << "Passo di campionamento: " << dt << " ns (" << 1.0 / dt << " GS/s)\n"
            << "Finestre [campioni]: amp " << AMP_LO_NS << "ns -" << AMP_HI_NS
            << "ns , baseline " << BN_LO_NS << "ns -" << BN_HI_NS
            << "ns , carica -" << Q_PRE_NS << "ns /+" << Q_POST_NS << "ns  attorno al minimo\n";

  TFile* fout = new TFile(outName, "RECREATE");
  TTree* tree = new TTree("events", "DRS4 events");

  UInt_t   b_id = 0;
  UShort_t b_tcell = 0;
  Double_t b_unixtime = 0;                 // secondi (fuso orario locale), con ms
  std::vector<double> amp(nch), charge(nch), bnoise(nch);

  tree->Branch("event_id", &b_id, "event_id/i");
  tree->Branch("trig_cell", &b_tcell, "trig_cell/s");
  tree->Branch("unixtime", &b_unixtime, "unixtime/D");
  for (size_t k = 0; k < nch; ++k) {
    int c = r.channels[k];
    tree->Branch(Form("ch%d_amp", c),    &amp[k],    Form("ch%d_amp/D", c));     // baseline - min [V]
    tree->Branch(Form("ch%d_charge", c), &charge[k], Form("ch%d_charge/D", c));  // [C], baseline sottratta
    tree->Branch(Form("ch%d_bnoise", c), &bnoise[k], Form("ch%d_bnoise/D", c));  // [V]
  }

  TH1F* hAmp = new TH1F("hAmp", Form("Pulse height, Channel %d;Ampiezza [mV];Eventi", ch),
                        NB_AMP, AMP_MIN_MV, AMP_MAX_MV);
  TH1F* hQ   = new TH1F("hQ", Form("Charge, Channel %d;Carica [pC];Eventi", ch),
                        NB_Q, Q_MIN_PC, Q_MAX_PC);

  DrsEvent ev;
  int n = 0;
  while (r.next(ev)) {
    for (size_t k = 0; k < nch; ++k) {
      int c = r.channels[k];
      std::vector<double> w(NCELL);
      for (int i = 0; i < NCELL; ++i) w[i] = ev.adc[c][i] / 65535.0 - 0.5;

      double m = 0, m2 = 0;                // baseline: media e std su [BN_LO, BN_HI)
      const int nb = BN_HI - BN_LO;
      for (int i = BN_LO; i < BN_HI; ++i) { m += w[i]; m2 += w[i] * w[i]; }
      m /= nb; m2 /= nb;
      bnoise[k] = std::sqrt(std::max(0.0, m2 - m * m));
      const double base = m;               // baseline = media su [BN_LO, BN_HI)

      // ampiezza: minimo in [AMP_LO, AMP_HI) rispetto alla baseline.
      // Il "-" serve per avere ampiezza positiva, visto che la waveform e' negativa -> eliminare in caso contrario
      const int imin = int(std::min_element(w.begin() + AMP_LO, w.begin() + AMP_HI) - w.begin());
      amp[k] = -(w[imin] - base);

      // carica: integrale di (w - baseline) in [imin - Q_PRE, imin + Q_POST)
      const int s0 = std::max(imin - Q_PRE, 0), s1 = std::min(imin + Q_POST, NCELL);
      double q = 0;
      for (int j = s0; j < s1; ++j)
        q += (w[j] - base) * r.tw.at(c)[(j - 1 + ev.trig_cell + NCELL) % NCELL] * 1e-9 / 50.0;
      charge[k] = -q;

      if (c == ch) { hAmp->Fill(amp[k] * 1e3); hQ->Fill(charge[k] * 1e12); }   // V -> mV, C -> pC
    }

    b_id = ev.id;
    b_tcell = ev.trig_cell;
    TDatime dtm(ev.date[0], ev.date[1], ev.date[2], ev.date[3], ev.date[4], ev.date[5]);
    b_unixtime = dtm.Convert() + ev.date[6] / 1000.0;
    tree->Fill();
    ++n;
  }
  std::cout << "Eventi letti: " << n << "\n";
  fout->Write();                           // scrive tree, hAmp e hQ

  if (hAmp->GetEntries() == 0) {
    std::cerr << "Canale " << ch << " non presente nel file (o nessun evento): niente plot.\n";
    fout->Close();
    return;
  }

  // gli istogrammi devono sopravvivere alla chiusura del file, perche' servono alle canvas
  hAmp->SetDirectory(nullptr);
  hQ->SetDirectory(nullptr);
  fout->Close();

  // stile (da impostare prima di creare le canvas)
  gStyle->SetOptStat("emr");               // entries, media, rms
  gStyle->SetPadGridX(true);  gStyle->SetPadGridY(true);
  gStyle->SetPadTickX(1);     gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.13);   gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.13); gStyle->SetPadTopMargin(0.09);
  gStyle->SetTitleFontSize(0.05);

  auto style = [](TH1F* h, int fill, int line) {
    h->SetFillColorAlpha(fill, 0.55);
    h->SetLineColor(line);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.045);  h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetTitleSize(0.045);  h->GetYaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetTitleOffset(1.1);  h->GetYaxis()->SetTitleOffset(1.3);
  };

  // canvas dell'ampiezza
  TCanvas* c1 = new TCanvas("c_amp", Form("Ampiezza, ch%d", ch), 900, 650);
  c1->SetLogy(LOGY);
  style(hAmp, kAzure + 1, kAzure + 3);
  hAmp->Draw("hist");
  c1->Update();
  c1->SaveAs(Form("amp_ch%d.png", ch));

  // canvas della carica
  TCanvas* c2 = new TCanvas("c_charge", Form("Carica, ch%d", ch), 900, 650);
  c2->SetLogy(LOGY);
  style(hQ, kOrange + 1, kOrange + 3);
  hQ->Draw("hist");
  c2->Update();
  c2->SaveAs(Form("charge_ch%d.png", ch));
}