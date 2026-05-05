// compareFieldMaps.C
//
// Compare the CERN-measured sPHENIX solenoid field map with the field map
// used in offline tracking analysis.
//
//   Measured:   sphenix-cernfinal-map/sPHENIXFieldMap  (CSV-based; mm, T)
//   Analysis:   fieldmap-used-in-analysis/...trackingmapxyz.root
//               TNtuple "fieldmap"; branches x,y,z,bx,by,bz,hz; units cm, T
//
// Run from the project directory:
//   root -l -b -q 'compareFieldMaps.C+'
//
// Major question: is the measured field consistent with the calculated
// solenoid field after a small rigid-body translation and/or rotation?
// Method: decompose (Bz_calc − Bz_meas), (Br_calc − Br_meas), and Bφ_calc
// into azimuthal Fourier modes.  A dominant m=1 signal with a phase that is
// CONSTANT across (r, z) is the signature of a rigid-body misalignment.

#include "sphenix-cernfinal-map/sPHENIXFieldMap.cxx"

#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TMath.h"
#include "TSystem.h"

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <string>

// ─── Analysis field map grid constants ───────────────────────────────────────
// Cartesian cube, x/y/z in cm, field in Tesla
static const int   NXC = 111, NYC = 111, NZC = 111;
static const float XMIN_C = -110.f, DX_C = 2.f;
static const float YMIN_C = -110.f, DY_C = 2.f;
static const float ZMIN_C = -110.f, DZ_C = 2.f;

// Heap-allocated 3-D grids (filled from the TNtuple)
static float *gCalcBx = nullptr;
static float *gCalcBy = nullptr;
static float *gCalcBz = nullptr;

inline int Cidx(int ix, int iy, int iz)
{ return ix * NYC*NZC + iy * NZC + iz; }

// ─── Trilinear interpolation of the analysis map ─────────────────────────────
// Returns false (and zeros) if the point is outside the grid.
bool CalcInterp(float xq, float yq, float zq,
                float &Bx, float &By, float &Bz)
{
    float fx = (xq - XMIN_C)/DX_C;
    float fy = (yq - YMIN_C)/DY_C;
    float fz = (zq - ZMIN_C)/DZ_C;
    if (fx < 0 || fx > NXC-1 || fy < 0 || fy > NYC-1 || fz < 0 || fz > NZC-1) {
        Bx = By = Bz = 0.f;
        return false;
    }
    // Clamp floor index so the right boundary (fx==NXC-1) maps to the last cell
    int ix = std::min((int)fx, NXC-2);
    int iy = std::min((int)fy, NYC-2);
    int iz = std::min((int)fz, NZC-2);
    float tx = fx-ix, ty = fy-iy, tz = fz-iz;
    auto I = [&](const float *G) -> float {
        return (G[Cidx(ix,  iy,  iz  )]*(1-tx)*(1-ty)*(1-tz)
              + G[Cidx(ix+1,iy,  iz  )]*   tx *(1-ty)*(1-tz)
              + G[Cidx(ix,  iy+1,iz  )]*(1-tx)*   ty *(1-tz)
              + G[Cidx(ix+1,iy+1,iz  )]*   tx *   ty *(1-tz)
              + G[Cidx(ix,  iy,  iz+1)]*(1-tx)*(1-ty)*   tz
              + G[Cidx(ix+1,iy,  iz+1)]*   tx *(1-ty)*   tz
              + G[Cidx(ix,  iy+1,iz+1)]*(1-tx)*   ty *   tz
              + G[Cidx(ix+1,iy+1,iz+1)]*   tx *   ty *   tz);
    };
    Bx = I(gCalcBx);  By = I(gCalcBy);  Bz = I(gCalcBz);
    return true;
}

// ─── Utility: save and delete a canvas ───────────────────────────────────────
void SaveClose(TCanvas *c, const char *fname)
{
    c->SaveAs(fname);
    delete c;
    printf("  -> %s\n", fname);
}

// ─── Utility: set a diverging blue-white-red palette ─────────────────────────
void UseBWR()
{
    const Int_t n = 5;
    Double_t st[n] = {0.00, 0.25, 0.50, 0.75, 1.00};
    Double_t r[n]  = {0.00, 0.00, 1.00, 1.00, 1.00};
    Double_t g[n]  = {0.00, 0.50, 1.00, 0.50, 0.00};
    Double_t b[n]  = {1.00, 1.00, 1.00, 0.00, 0.00};
    TColor::CreateGradientColorTable(n, st, r, g, b, 255);
    gStyle->SetNumberContours(255);
}

// ─── Utility: set symmetric Z range for a diverging 2-D histogram ────────────
void SymZ(TH2F *h, double lim = -1.)
{
    if (lim < 0) lim = std::max(std::abs(h->GetMaximum()),
                                std::abs(h->GetMinimum()));
    h->GetZaxis()->SetRangeUser(-lim, lim);
}


// ═════════════════════════════════════════════════════════════════════════════
void compareFieldMaps()
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetPadRightMargin(0.13);
    gSystem->mkdir("plots", kTRUE);

    // ── Paths ─────────────────────────────────────────────────────────────────
    const char *FINE_CSV  =
        "/Users/haggerty/sphenix/data/data02/sphenix/MagnetMapping/cernfinal"
        "/fieldMapFineFullField.csv";
    const char *ROUGH_CSV =
        "/Users/haggerty/sphenix/data/data02/sphenix/MagnetMapping/cernfinal"
        "/fieldMapRoughFullField.csv";
    const char *CALC_ROOT =
        "fieldmap-used-in-analysis/"
        "8e4d6c3b1660540a658da3a275af2bde_sphenix3dtrackingmapxyz.root";

    // ─────────────────────────────────────────────────────────────────────────
    // 1. Load measured field map
    // ─────────────────────────────────────────────────────────────────────────
    printf("=== Loading measured field map ===\n");
    sPHENIXFieldMap meas(FINE_CSV, ROUGH_CSV);

    // Grid constants (from sPHENIXFieldMap.h)
    const int    NR_M = 37, NZ_M = 241;
    const double DR_MM = 25., DZ_MM = 20., ZMIN_MM = -2700.;

    double bz_m[NR_M][NZ_M], br_m[NR_M][NZ_M];
    for (int ir = 0; ir < NR_M; ++ir)
        for (int iz = 0; iz < NZ_M; ++iz) {
            bz_m[ir][iz] = meas.GetBzGrid(ir, iz);
            br_m[ir][iz] = meas.GetBrGrid(ir, iz);
        }
    printf("  Read %d x %d measured-map grid.\n", NR_M, NZ_M);

    // ─────────────────────────────────────────────────────────────────────────
    // 2. Load analysis field map into 3-D Cartesian arrays
    // ─────────────────────────────────────────────────────────────────────────
    printf("=== Loading analysis field map ===\n");
    TFile *fCalc = TFile::Open(CALC_ROOT, "READ");
    if (!fCalc || fCalc->IsZombie()) {
        fprintf(stderr, "ERROR: cannot open %s\n", CALC_ROOT); return;
    }
    TNtuple *nt = (TNtuple*)fCalc->Get("fieldmap");
    Long64_t nEnt = nt->GetEntries();
    printf("  %lld entries\n", nEnt);

    gCalcBx = new float[NXC*NYC*NZC]();
    gCalcBy = new float[NXC*NYC*NZC]();
    gCalcBz = new float[NXC*NYC*NZC]();

    float xb, yb, zb, bxb, byb, bzb;
    nt->SetBranchAddress("x",  &xb);
    nt->SetBranchAddress("y",  &yb);
    nt->SetBranchAddress("z",  &zb);
    nt->SetBranchAddress("bx", &bxb);
    nt->SetBranchAddress("by", &byb);
    nt->SetBranchAddress("bz", &bzb);

    for (Long64_t i = 0; i < nEnt; ++i) {
        nt->GetEntry(i);
        int ix = (int)std::lround((xb - XMIN_C)/DX_C);
        int iy = (int)std::lround((yb - YMIN_C)/DY_C);
        int iz = (int)std::lround((zb - ZMIN_C)/DZ_C);
        if (ix < 0 || ix >= NXC || iy < 0 || iy >= NYC || iz < 0 || iz >= NZC) continue;
        gCalcBx[Cidx(ix,iy,iz)] = bxb;
        gCalcBy[Cidx(ix,iy,iz)] = byb;
        gCalcBz[Cidx(ix,iy,iz)] = bzb;
    }
    fCalc->Close();
    printf("  3-D grid filled.\n");

    // ─────────────────────────────────────────────────────────────────────────
    // 3. Comparison grid
    //
    //   r:   0, 2.5, 5, ..., 75 cm  (NR=31, matches measured-map r-step)
    //   phi: 36 azimuths, 10-deg apart
    //   z:   −110 to +110 cm, 2-cm steps  (NZ=111)
    //        → overlap of measured (−270..+210 cm) and analysis (−110..+110 cm)
    //        → 2-cm step matches both grids exactly
    // ─────────────────────────────────────────────────────────────────────────
    const int    NR  = 31;    // r = 0, 2.5, ..., 75 cm
    const int    NZ  = 111;   // z = −110, −108, ..., +110 cm
    const int    NP  = 36;    // phi = 0, 10, ..., 350 deg
    const double DR  = 2.5;   // cm
    const double DZ  = 2.0;   // cm
    const double Z0  = -110.; // cm, start of overlap

    // Index into measured-map z-array for the overlap start
    const int IZ0 = (int)std::lround((Z0*10. - ZMIN_MM) / DZ_MM);  // = 80
    printf("  Overlap: measured iz = %d .. %d  (z = %.0f .. %.0f mm)\n",
           IZ0, IZ0+NZ-1, ZMIN_MM + IZ0*DZ_MM, ZMIN_MM + (IZ0+NZ-1)*DZ_MM);

    // ─────────────────────────────────────────────────────────────────────────
    // 4. Evaluate analysis map on (r, phi, z) grid and accumulate Fourier modes
    //
    //   m=0:  Fx_m0[ir][iz] = (1/NP) Σ_phi F(r, phi, z)
    //   m=1:  Fx_a1[ir][iz] = (2/NP) Σ_phi F(r, phi, z) cos(phi)
    //         Fx_b1[ir][iz] = (2/NP) Σ_phi F(r, phi, z) sin(phi)
    //   A1   = sqrt(a1^2 + b1^2)
    //   For a tilted/offset solenoid the m=1 phase phi1 = atan2(b1,a1)
    //   should be the SAME everywhere in (r, z).
    // ─────────────────────────────────────────────────────────────────────────
    printf("=== Evaluating analysis map on %d r x %d phi x %d z points ===\n",
           NR, NP, NZ);

    // Result 2-D arrays [ir][iz]
    double Bz_m0[NR][NZ]  = {}, Bz_a1[NR][NZ]  = {}, Bz_b1[NR][NZ]  = {};
    double Br_m0[NR][NZ]  = {}, Br_a1[NR][NZ]  = {}, Br_b1[NR][NZ]  = {};
    double Bp_m0[NR][NZ]  = {}, Bp_a1[NR][NZ]  = {}, Bp_b1[NR][NZ]  = {};

    for (int ir = 0; ir < NR; ++ir) {
        float r = (float)(ir * DR);
        for (int ip = 0; ip < NP; ++ip) {
            double phi = ip * (TMath::TwoPi() / NP);
            float  xq  = r * (float)std::cos(phi);
            float  yq  = r * (float)std::sin(phi);
            double cp  = std::cos(phi), sp = std::sin(phi);
            for (int iz = 0; iz < NZ; ++iz) {
                float zq = (float)(Z0 + iz*DZ);
                float Bx, By, Bz;
                CalcInterp(xq, yq, zq, Bx, By, Bz);
                double Br  =  Bx*cp + By*sp;   // cylindrical radial
                double Bphi= -Bx*sp + By*cp;   // cylindrical azimuthal
                Bz_m0[ir][iz] += Bz;   Bz_a1[ir][iz] += Bz*cp;  Bz_b1[ir][iz] += Bz*sp;
                Br_m0[ir][iz] += Br;   Br_a1[ir][iz] += Br*cp;  Br_b1[ir][iz] += Br*sp;
                Bp_m0[ir][iz] += Bphi; Bp_a1[ir][iz] += Bphi*cp;Bp_b1[ir][iz] += Bphi*sp;
            }
        }
    }
    // Normalize
    for (int ir = 0; ir < NR; ++ir)
        for (int iz = 0; iz < NZ; ++iz) {
            Bz_m0[ir][iz] /= NP;  Bz_a1[ir][iz] *= 2./NP;  Bz_b1[ir][iz] *= 2./NP;
            Br_m0[ir][iz] /= NP;  Br_a1[ir][iz] *= 2./NP;  Br_b1[ir][iz] *= 2./NP;
            Bp_m0[ir][iz] /= NP;  Bp_a1[ir][iz] *= 2./NP;  Bp_b1[ir][iz] *= 2./NP;
        }

    // ─────────────────────────────────────────────────────────────────────────
    // 5. Differences and amplitudes
    // ─────────────────────────────────────────────────────────────────────────
    double dBz[NR][NZ], dBr[NR][NZ];
    double Bz_A1[NR][NZ], Br_A1[NR][NZ], Bp_A1[NR][NZ];
    double Bz_ph[NR][NZ], Bp_ph[NR][NZ];

    double maxdBz=0., maxdBr=0., maxBp_m0=0.;
    double maxBz_A1=0., maxBr_A1=0., maxBp_A1=0.;

    for (int ir = 0; ir < NR; ++ir)
        for (int iz = 0; iz < NZ; ++iz) {
            dBz[ir][iz] = Bz_m0[ir][iz] - bz_m[ir][IZ0+iz];
            dBr[ir][iz] = Br_m0[ir][iz] - br_m[ir][IZ0+iz];
            Bz_A1[ir][iz] = std::hypot(Bz_a1[ir][iz], Bz_b1[ir][iz]);
            Br_A1[ir][iz] = std::hypot(Br_a1[ir][iz], Br_b1[ir][iz]);
            Bp_A1[ir][iz] = std::hypot(Bp_a1[ir][iz], Bp_b1[ir][iz]);
            Bz_ph[ir][iz] = std::atan2(Bz_b1[ir][iz], Bz_a1[ir][iz]) * TMath::RadToDeg();
            Bp_ph[ir][iz] = std::atan2(Bp_b1[ir][iz], Bp_a1[ir][iz]) * TMath::RadToDeg();
            maxdBz  = std::max(maxdBz,  std::abs(dBz[ir][iz]));
            maxdBr  = std::max(maxdBr,  std::abs(dBr[ir][iz]));
            maxBp_m0= std::max(maxBp_m0,std::abs(Bp_m0[ir][iz]));
            maxBz_A1= std::max(maxBz_A1,Bz_A1[ir][iz]);
            maxBr_A1= std::max(maxBr_A1,Br_A1[ir][iz]);
            maxBp_A1= std::max(maxBp_A1,Bp_A1[ir][iz]);
        }

    // ─────────────────────────────────────────────────────────────────────────
    // 6. Fill TH2F histograms  (x-axis = z, y-axis = r)
    // ─────────────────────────────────────────────────────────────────────────
    // Histogram axes: z in cm (x), r in cm (y)
    // Bin edges: z from Z0 to Z0+(NZ)*DZ, r from 0 to NR*DR
    const double ZMAX = Z0 + NZ*DZ;
    const double RMAX = NR*DR;

    auto MH2 = [&](const char *nm, const char *ttl) -> TH2F* {
        return new TH2F(nm, ttl, NZ, Z0, ZMAX, NR, 0., RMAX);
    };

    TH2F *hBz_meas = MH2("hBz_meas", "Measured Bz;z (cm);r (cm)");
    TH2F *hBz_calc = MH2("hBz_calc", "Analysis Bz  #phi-avg;z (cm);r (cm)");
    TH2F *hdBz     = MH2("hdBz",     "#DeltaBz (Analysis #minus Measured) [mT];z (cm);r (cm)");
    TH2F *hBr_meas = MH2("hBr_meas", "Measured Br;z (cm);r (cm)");
    TH2F *hBr_calc = MH2("hBr_calc", "Analysis Br  #phi-avg;z (cm);r (cm)");
    TH2F *hdBr     = MH2("hdBr",     "#DeltaBr (Analysis #minus Measured) [mT];z (cm);r (cm)");
    TH2F *hBp_m0   = MH2("hBp_m0",  "B#phi  #phi-mean [mT];z (cm);r (cm)");
    TH2F *hBp_A1   = MH2("hBp_A1",  "B#phi  m=1 amplitude [mT];z (cm);r (cm)");
    TH2F *hBz_A1h  = MH2("hBz_A1",  "Bz  m=1 amplitude [mT];z (cm);r (cm)");
    TH2F *hBr_A1h  = MH2("hBr_A1",  "Br  m=1 amplitude [mT];z (cm);r (cm)");
    TH2F *hBp_ph   = MH2("hBp_ph",  "B#phi  m=1 phase [deg];z (cm);r (cm)");
    TH2F *hBz_ph_h = MH2("hBz_ph",  "Bz  m=1 phase [deg];z (cm);r (cm)");

    for (int ir = 0; ir < NR; ++ir)
        for (int iz = 0; iz < NZ; ++iz) {
            // SetBinContent(ix, iy, val): ix = z-bin (1-indexed), iy = r-bin
            int bx = iz+1, by = ir+1;
            hBz_meas->SetBinContent(bx, by, bz_m[ir][IZ0+iz]);
            hBz_calc->SetBinContent(bx, by, Bz_m0[ir][iz]);
            hdBz    ->SetBinContent(bx, by, dBz[ir][iz]*1e3);
            hBr_meas->SetBinContent(bx, by, br_m[ir][IZ0+iz]);
            hBr_calc->SetBinContent(bx, by, Br_m0[ir][iz]);
            hdBr    ->SetBinContent(bx, by, dBr[ir][iz]*1e3);
            hBp_m0  ->SetBinContent(bx, by, Bp_m0[ir][iz]*1e3);
            hBp_A1  ->SetBinContent(bx, by, Bp_A1[ir][iz]*1e3);
            hBz_A1h ->SetBinContent(bx, by, Bz_A1[ir][iz]*1e3);
            hBr_A1h ->SetBinContent(bx, by, Br_A1[ir][iz]*1e3);
            // Phase: set to 0 in low-amplitude cells (will be visually ambiguous
            // in those regions, but amplitude plots are shown alongside)
            hBp_ph  ->SetBinContent(bx, by, Bp_ph[ir][iz]);
            hBz_ph_h->SetBinContent(bx, by, Bz_ph[ir][iz]);
        }

    // ─────────────────────────────────────────────────────────────────────────
    // 7. Plots
    // ─────────────────────────────────────────────────────────────────────────
    printf("=== Making plots ===\n");

    // ── 01. On-axis Bz(z) comparison ─────────────────────────────────────────
    {
        TGraph *gM  = new TGraph(NZ);
        TGraph *gC  = new TGraph(NZ);
        TGraph *gD  = new TGraph(NZ);
        for (int iz = 0; iz < NZ; ++iz) {
            double z = Z0 + iz*DZ;
            gM->SetPoint(iz, z, bz_m[0][IZ0+iz]);
            gC->SetPoint(iz, z, Bz_m0[0][iz]);
            gD->SetPoint(iz, z, dBz[0][iz]*1e3);
        }
        gM->SetLineColor(kBlue);  gM->SetLineWidth(2);
        gC->SetLineColor(kRed);   gC->SetLineWidth(2); gC->SetLineStyle(2);
        gD->SetLineColor(kBlack); gD->SetLineWidth(2);

        TCanvas *c = new TCanvas("c01", "", 900, 700);
        c->Divide(1, 2);

        c->cd(1);
        TMultiGraph *mg = new TMultiGraph("mg_onaxis",
            "On-axis Bz(r=0, z);z (cm);Bz (T)");
        mg->Add(gM, "L"); mg->Add(gC, "L");
        mg->Draw("A");
        TLegend *leg = new TLegend(0.65, 0.15, 0.88, 0.35);
        leg->AddEntry(gM, "Measured", "l");
        leg->AddEntry(gC, "Analysis (#phi-avg)", "l");
        leg->Draw();

        c->cd(2);
        gD->SetTitle("On-axis #DeltaBz = Analysis #minus Measured;z (cm);#DeltaBz (mT)");
        gD->Draw("AL");
        TLine *zl = new TLine(Z0, 0., ZMAX, 0.);
        zl->SetLineStyle(2); zl->SetLineColor(kGray+1); zl->Draw();

        SaveClose(c, "plots/01_onaxis_Bz.pdf");
    }

    // ── 02. 2-D Bz maps ───────────────────────────────────────────────────────
    {
        gStyle->SetPalette(kBird);
        TCanvas *c = new TCanvas("c02", "", 1400, 500);
        c->Divide(2, 1);
        hBz_meas->GetZaxis()->SetRangeUser(-0.1, 1.65);
        hBz_calc->GetZaxis()->SetRangeUser(-0.1, 1.65);
        c->cd(1); hBz_meas->Draw("COLZ");
        c->cd(2); hBz_calc->Draw("COLZ");
        SaveClose(c, "plots/02_Bz_2d_maps.pdf");
    }

    // ── 03. ΔBz m=0 ───────────────────────────────────────────────────────────
    {
        UseBWR();
        SymZ(hdBz);
        TCanvas *c = new TCanvas("c03", "", 950, 500);
        hdBz->Draw("COLZ");
        SaveClose(c, "plots/03_dBz_m0.pdf");
    }

    // ── 04. Br comparison ─────────────────────────────────────────────────────
    {
        UseBWR();
        SymZ(hBr_meas); SymZ(hBr_calc); SymZ(hdBr);
        TCanvas *c = new TCanvas("c04", "", 1700, 500);
        c->Divide(3, 1);
        c->cd(1); hBr_meas->Draw("COLZ");
        c->cd(2); hBr_calc->Draw("COLZ");
        c->cd(3); hdBr->Draw("COLZ");
        SaveClose(c, "plots/04_Br_comparison.pdf");
    }

    // ── 05. m=1 amplitude maps ────────────────────────────────────────────────
    // Non-zero amplitudes at consistent phi → rigid-body solenoid misalignment
    {
        gStyle->SetPalette(kViridis);
        TCanvas *c = new TCanvas("c05", "", 1700, 500);
        c->Divide(3, 1);
        c->cd(1); hBz_A1h->Draw("COLZ");
        c->cd(2); hBr_A1h->Draw("COLZ");
        c->cd(3); hBp_A1->Draw("COLZ");
        c->SetTitle("m=1 amplitudes: constant phase across (r,z) "
                    "=> rigid-body misalignment");
        SaveClose(c, "plots/05_m1_amplitudes.pdf");
    }

    // ── 06. m=1 phase maps (amplitude + phase side by side) ───────────────────
    // A phase that is CONSTANT over (r,z) indicates a single dipole direction.
    {
        TCanvas *c = new TCanvas("c06", "", 1400, 900);
        c->Divide(2, 2);

        gStyle->SetPalette(kViridis);
        c->cd(1); hBp_A1->Draw("COLZ");

        gStyle->SetPalette(kRainBow);
        hBp_ph->GetZaxis()->SetRangeUser(-180., 180.);
        c->cd(2); hBp_ph->Draw("COLZ");

        gStyle->SetPalette(kViridis);
        c->cd(3); hBz_A1h->Draw("COLZ");

        gStyle->SetPalette(kRainBow);
        hBz_ph_h->GetZaxis()->SetRangeUser(-180., 180.);
        c->cd(4); hBz_ph_h->Draw("COLZ");

        SaveClose(c, "plots/06_m1_amplitude_and_phase.pdf");
    }

    // ── 07. Bphi overview (m=0 mean and m=1 amplitude) ────────────────────────
    {
        UseBWR();
        SymZ(hBp_m0);
        TCanvas *c = new TCanvas("c07", "", 1400, 500);
        c->Divide(2, 1);
        c->cd(1); hBp_m0->Draw("COLZ");
        gStyle->SetPalette(kViridis);
        c->cd(2); hBp_A1->Draw("COLZ");
        SaveClose(c, "plots/07_Bphi_overview.pdf");
    }

    // ── 08. Bz radial profiles at selected z values ───────────────────────────
    {
        const int NSZ = 5;
        double zsamp[NSZ] = {0., 50., 100., -50., -100.};
        TCanvas *c = new TCanvas("c08", "", 1800, 420);
        c->Divide(NSZ, 1);
        for (int s = 0; s < NSZ; ++s) {
            c->cd(s+1);
            int iz_s = (int)std::lround((zsamp[s] - Z0)/DZ);
            TGraph *gMr = new TGraph(NR), *gCr = new TGraph(NR);
            for (int ir = 0; ir < NR; ++ir) {
                gMr->SetPoint(ir, ir*DR, bz_m[ir][IZ0+iz_s]);
                gCr->SetPoint(ir, ir*DR, Bz_m0[ir][iz_s]);
            }
            gMr->SetLineColor(kBlue); gMr->SetLineWidth(2);
            gCr->SetLineColor(kRed);  gCr->SetLineWidth(2); gCr->SetLineStyle(2);
            TString ttl = Form("z = %+.0f cm;r (cm);Bz (T)", zsamp[s]);
            TMultiGraph *mg = new TMultiGraph(Form("mg_r%d", s), ttl);
            mg->Add(gMr, "L"); mg->Add(gCr, "L");
            mg->Draw("A");
            if (s == 0) {
                TLegend *leg = new TLegend(0.45, 0.7, 0.98, 0.92);
                leg->AddEntry(gMr, "Measured", "l");
                leg->AddEntry(gCr, "Analysis #phi-avg", "l");
                leg->Draw();
            }
        }
        SaveClose(c, "plots/08_Bz_radial_profiles.pdf");
    }

    // ── 09. Field vs phi at a sample tracking point (r=30 cm, z=0) ───────────
    // Sinusoidal variation at m=1 → rigid-body misalignment
    // Higher harmonics → non-trivial field distortion
    {
        const double R_SAMP = 30., Z_SAMP = 0.;
        int ir_s = (int)std::lround(R_SAMP/DR);
        int iz_s = (int)std::lround((Z_SAMP - Z0)/DZ);

        TGraph *gBzP = new TGraph(NP);
        TGraph *gBrP = new TGraph(NP);
        TGraph *gBpP = new TGraph(NP);
        for (int ip = 0; ip < NP; ++ip) {
            double phi = ip * (TMath::TwoPi()/NP);
            double cp = std::cos(phi), sp = std::sin(phi);
            float xq = (float)(R_SAMP*cp), yq = (float)(R_SAMP*sp);
            float zq = (float)(Z0 + iz_s*DZ);
            float Bx, By, Bz;
            CalcInterp(xq, yq, zq, Bx, By, Bz);
            gBzP->SetPoint(ip, ip*360./NP, Bz);
            gBrP->SetPoint(ip, ip*360./NP,  Bx*cp + By*sp);
            gBpP->SetPoint(ip, ip*360./NP, -Bx*sp + By*cp);
        }
        gBzP->SetLineColor(kRed); gBzP->SetLineWidth(2);
        gBrP->SetLineColor(kRed); gBrP->SetLineWidth(2);
        gBpP->SetLineColor(kRed); gBpP->SetLineWidth(2);

        TCanvas *c = new TCanvas("c09", "", 1500, 430);
        c->Divide(3, 1);

        auto DrawPhi = [&](int pad, TGraph *g, double measVal,
                           const char *yl) {
            c->cd(pad);
            g->SetTitle(Form("%s vs #phi  (r=%.0f cm, z=%.0f cm);#phi (deg);%s",
                             yl, R_SAMP, Z_SAMP, yl));
            g->Draw("AL");
            TLine *ml = new TLine(0., measVal, 350., measVal);
            ml->SetLineColor(kBlue); ml->SetLineWidth(2); ml->SetLineStyle(2);
            ml->Draw();
            TLegend *leg = new TLegend(0.58, 0.77, 0.98, 0.92);
            leg->AddEntry(g,  "Analysis",            "l");
            leg->AddEntry(ml, "Measured (#phi-avg)", "l");
            leg->Draw();
        };
        DrawPhi(1, gBzP, bz_m[ir_s][IZ0+iz_s], "Bz (T)");
        DrawPhi(2, gBrP, br_m[ir_s][IZ0+iz_s], "Br (T)");
        DrawPhi(3, gBpP, 0.,                    "B#phi (T)");
        SaveClose(c, "plots/09_phi_dependence_r30_z0.pdf");
    }

    // ── 10. ΔBz at phi=0 versus phi=180 deg ───────────────────────────────────
    // For a transverse offset of the solenoid axis, ΔBz will be equal and
    // opposite at phi=0 and phi=180.  This plot makes that anti-symmetry visible.
    {
        TH2F *hdBz0   = MH2("hdBz0",   "#DeltaBz at #phi=0#circ [mT];z (cm);r (cm)");
        TH2F *hdBz180 = MH2("hdBz180", "#DeltaBz at #phi=180#circ [mT];z (cm);r (cm)");
        for (int ir = 0; ir < NR; ++ir)
            for (int iz = 0; iz < NZ; ++iz) {
                float zq = (float)(Z0 + iz*DZ);
                float Bx0, By0, Bz0, Bx180, By180, Bz180;
                CalcInterp((float)(ir*DR), 0.f, zq, Bx0, By0, Bz0);
                CalcInterp((float)(-ir*DR), 0.f, zq, Bx180, By180, Bz180);
                hdBz0  ->SetBinContent(iz+1, ir+1, (Bz0   - bz_m[ir][IZ0+iz])*1e3);
                hdBz180->SetBinContent(iz+1, ir+1, (Bz180 - bz_m[ir][IZ0+iz])*1e3);
            }
        UseBWR();
        SymZ(hdBz0); SymZ(hdBz180);
        TCanvas *c = new TCanvas("c10", "", 1400, 500);
        c->Divide(2, 1);
        c->cd(1); hdBz0->Draw("COLZ");
        c->cd(2); hdBz180->Draw("COLZ");
        SaveClose(c, "plots/10_dBz_phi0_vs_phi180.pdf");
        delete hdBz0; delete hdBz180;
    }

    // ─────────────────────────────────────────────────────────────────────────
    // 8. Save all histograms to a ROOT file
    // ─────────────────────────────────────────────────────────────────────────
    TFile *fOut = TFile::Open("plots/comparison_histograms.root", "RECREATE");
    hBz_meas->Write(); hBz_calc->Write(); hdBz->Write();
    hBr_meas->Write(); hBr_calc->Write(); hdBr->Write();
    hBp_m0->Write();   hBp_A1->Write();   hBp_ph->Write();
    hBz_A1h->Write();  hBz_ph_h->Write();
    hBr_A1h->Write();
    fOut->Close();
    printf("  -> plots/comparison_histograms.root\n");

    // ─────────────────────────────────────────────────────────────────────────
    // 9. Numerical summary
    // ─────────────────────────────────────────────────────────────────────────
    int iz_z0 = (int)std::lround((0. - Z0)/DZ);
    printf("\n=== Summary ===\n");
    printf("  On-axis Bz(z=0):   measured = %.4f T,  analysis = %.4f T,  "
           "diff = %+.1f mT\n",
           bz_m[0][IZ0+iz_z0], Bz_m0[0][iz_z0], dBz[0][iz_z0]*1e3);
    printf("  Max |ΔBz| m=0  in tracking volume : %.1f mT\n", maxdBz*1e3);
    printf("  Max |ΔBr| m=0  in tracking volume : %.1f mT\n", maxdBr*1e3);
    printf("  Max |Bφ| azimuthal mean            : %.2f mT\n", maxBp_m0*1e3);
    printf("  Max |Bz| m=1 amplitude             : %.2f mT\n", maxBz_A1*1e3);
    printf("  Max |Br| m=1 amplitude             : %.2f mT\n", maxBr_A1*1e3);
    printf("  Max |Bφ| m=1 amplitude             : %.2f mT\n", maxBp_A1*1e3);

    // Circular mean and spread of the Bphi m=1 phase (significant cells only)
    {
        double sx = 0., sy = 0., n = 0.;
        double thresh = 0.3 * maxBp_A1;
        for (int ir = 0; ir < NR; ++ir)
            for (int iz = 0; iz < NZ; ++iz)
                if (Bp_A1[ir][iz] > thresh) {
                    sx += std::cos(Bp_ph[ir][iz] * TMath::DegToRad());
                    sy += std::sin(Bp_ph[ir][iz] * TMath::DegToRad());
                    n += 1.;
                }
        if (n > 0.) {
            double mean_ph = std::atan2(sy/n, sx/n) * TMath::RadToDeg();
            double R_circ  = std::hypot(sx/n, sy/n);
            double spread  = (R_circ > 0.) ?
                             std::sqrt(-2.*std::log(R_circ)) * TMath::RadToDeg() : 999.;
            printf("  Bφ m=1 phase (where A1 > 30%% of max): "
                   "%.1f +/- %.1f deg  (%.0f cells)\n",
                   mean_ph, spread, n);
            printf("  Small spread => consistent with single tilt/offset direction\n");
        }
    }

    // Clean up
    delete[] gCalcBx; delete[] gCalcBy; delete[] gCalcBz;
    gCalcBx = gCalcBy = gCalcBz = nullptr;
    printf("\nDone.\n");
}
