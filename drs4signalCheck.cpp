// drs4_check.C
// Disegna eventi di un file binario DRS4, con tutti i canali sovrapposti e le
// finestre di analisi (baseline in grigio, ampiezza in rosso).
//  - modalita' normale: N eventi casuali in un'unica canvas (un pad per evento)
//  - modalita' opzionale: eventi scelti per numero (event_id), una canvas per evento
//
// Uso (senza -b e senza -q, altrimenti le finestre non restano aperte):
//   root -l 'drs4_check.C+("run.dat", 9, 0)'
//   root -l 'drs4_check.C+("run.dat", 9, 0, -0.5, 0.5, "12,45,300")'
//   (argomenti: file binario, n. eventi casuali, seed [0 = casuale], ymin, ymax,
//    numeri di evento separati da virgola [opzionale: se presenti si disegnano solo questi])

#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "TROOT.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TRandom3.h"

constexpr int NCELL = 1024;
// stesse finestre di drs4_plot.C (IN NANOSECONDI): modificale qui se le cambi la'
constexpr double AMP_LO_NS = 56.0, AMP_HI_NS = 70.0;   // ampiezza (-min)
constexpr double BN_LO_NS  = 2.0,  BN_HI_NS  = 54.0;   // baseline
double ymin = -0.5; double ymax = 0.5;   // override dei parametri di input, perche' i canali sono in V

struct DrsEvent {
  uint32_t id = 0;
  uint16_t date[7] = {0};
  uint16_t range_center = 0;
  uint16_t trig_cell = 0;
  std::map<int, uint32_t> scaler;
  std::map<int, std::vector<uint16_t>> adc;
};

class DrsReader {
public:
  std::vector<int> channels;
  std::map<int, std::vector<float>> tw;
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
    f.seekg(-4, std::ios::cur);
    return bool(f);
  }

  bool next(DrsEvent& ev) {
    char h[4];
    if (!f.read(h, 4)) return false;
    if (std::strncmp(h, "EHDR", 4)) return false;
    f.read(reinterpret_cast<char*>(&ev.id), 4);
    f.read(reinterpret_cast<char*>(ev.date), 14);
    f.read(reinterpret_cast<char*>(&ev.range_center), 2);
    f.read(h, 2);
    uint16_t bid; f.read(reinterpret_cast<char*>(&bid), 2);
    f.read(h, 2);
    f.read(reinterpret_cast<char*>(&ev.trig_cell), 2);
    for (int ch : channels) {
      f.read(h, 4);
      f.read(reinterpret_cast<char*>(&ev.scaler[ch]), 4);
      ev.adc[ch].resize(NCELL);
      f.read(reinterpret_cast<char*>(ev.adc[ch].data()), NCELL * 2);
    }
    return bool(f);
  }

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

void drs4_check(const char* fname = "run.dat", int nPick = 9, unsigned seed = 0, const char* evts = "") {
  //   gROOT->SetBatch(kTRUE);        // Questo non fa popuppare il plot

  DrsReader r;
  if (!r.open(fname)) return;
  if (r.channels.empty()) return;

  // passo di campionamento medio [ns] dalle time width, poi ns -> indice di campione
  double sumTw = 0;
  for (int c : r.channels) for (float w : r.tw[c]) sumTw += w;
  const double dt = sumTw / (double(NCELL) * r.channels.size());
  auto ns2s = [&](double ns) { return std::min(NCELL, std::max(0, (int)std::lround(ns / dt))); };
  const int AMP_LO = ns2s(AMP_LO_NS), AMP_HI = ns2s(AMP_HI_NS);
  const int BN_LO  = ns2s(BN_LO_NS),  BN_HI  = ns2s(BN_HI_NS);

  // numeri di evento richiesti (opzionale), es. "12,45,300"
  std::set<uint32_t> wanted;
  {
    std::stringstream ss(evts);
    std::string tok;
    while (std::getline(ss, tok, ','))
      if (!tok.empty()) wanted.insert((uint32_t)std::stoul(tok));
  }

  const int colors[4] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};

  // disegna un evento nel pad corrente: tutti i canali sovrapposti, legenda, finestre di analisi
  auto draw_event = [&](const DrsEvent& e) {
    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = new TLegend(0.78, 0.70, 0.98, 0.90);
    std::vector<double> tref;

    for (size_t k = 0; k < r.channels.size(); ++k) {
      int ch = r.channels[k];
      std::vector<double> t = r.time(ch, e.trig_cell);
      std::vector<double> v(NCELL);
      for (int j = 0; j < NCELL; ++j) v[j] = e.adc.at(ch)[j] / 65535.0 - 0.5;
      if (tref.empty()) tref = t;

      TGraph* g = new TGraph(NCELL, t.data(), v.data());
      g->SetLineColor(colors[k % 4]);
      g->SetLineWidth(1);
      mg->Add(g, "L");
      leg->AddEntry(g, Form("ch%d", ch), "l");
    }
    
    mg->SetTitle(Form("evento %u, trig cell %u;Time [ns];Voltage [V]", e.id, e.trig_cell));
    mg->SetMinimum(ymin);
    mg->SetMaximum(ymax);
    mg->Draw("AL");
    leg->Draw();

    // finestre di analisi (sull'asse tempi del primo canale)
    auto vline = [&](int cell, int col) {
      TLine* l = new TLine(tref[cell], ymin, tref[cell], ymax);
      l->SetLineColor(col);
      l->SetLineStyle(2);
      l->Draw();
    };
    vline(BN_LO, kGray + 1);   vline(BN_HI, kGray + 1);
    vline(AMP_LO, kRed);       vline(AMP_HI, kRed);
  };

  DrsEvent ev;

  // ---- modalita' opzionale: eventi scelti per numero, una canvas per evento ----
  if (!wanted.empty()) {
    std::vector<DrsEvent> found;
    std::set<uint32_t> foundIds;
    while (r.next(ev)) {
      if (wanted.count(ev.id) && foundIds.insert(ev.id).second) found.push_back(ev);
      if (foundIds.size() == wanted.size()) break;       // trovati tutti: smetto di leggere
    }
    for (uint32_t id : wanted)
      if (!foundIds.count(id)) std::cerr << "Evento " << id << " non trovato nel file.\n";
    if (found.empty()) return;

    for (const DrsEvent& e : found) {
      TCanvas* cv = new TCanvas(Form("ev%u", e.id), Form("Evento %u", e.id), 1200, 700);
      cv->cd();
      draw_event(e);
      cv->Update();
    }
    return;
  }

  // ---- modalita' normale: eventi casuali in un'unica canvas ----
  // reservoir sampling: un solo passaggio sul file
  TRandom3 rng(seed);                       // seed 0 -> casuale a ogni lancio
  std::vector<DrsEvent> pick;
  pick.reserve(nPick);
  long n = 0;
  while (r.next(ev)) {
    if ((int)pick.size() < nPick) pick.push_back(ev);
    else {
      long j = rng.Integer(n + 1);
      if (j < nPick) pick[j] = ev;
    }
    ++n;
  }
  std::cout << "Eventi nel file: " << n << ", disegnati: " << pick.size() << "\n";
  if (pick.empty()) return;
  std::sort(pick.begin(), pick.end(),
            [](const DrsEvent& a, const DrsEvent& b) { return a.id < b.id; });

  const int ncols = 3;
  const int nrows = ((int)pick.size() + ncols - 1) / ncols;
  TCanvas* c = new TCanvas("c", "DRS4 check", 1500, 300 * nrows);
  c->Divide(ncols, nrows);

  for (size_t i = 0; i < pick.size(); ++i) {
    c->cd(i + 1);
    // gPad->SetGrid();
    draw_event(pick[i]);
  }

  //   c->SaveAs("check_events.png");
}