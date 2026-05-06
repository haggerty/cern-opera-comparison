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
// Major question: is the measured field consistent with the OPERA
// solenoid field after a small rigid-body translation and/or rotation?
// Method: decompose (Bz_OPERA − Bz_meas), (Br_OPERA − Br_meas), and Bφ_OPERA
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

// ─── Trilinear interpolation of the OPERA map ────────────────────────────────
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
// calcRoot  – path to the Cartesian TNtuple to compare against the measured map.
//             Default: the OPERA file in fieldmap-used-in-analysis/.
// calcLabel – short name used in histogram titles and legend entries.
//             Default: "OPERA"
// outDir    – directory for output PDFs and ROOT file.
//             Default: "plots"
void compareFieldMaps(const char *calcRoot  = nullptr,
                      const char *calcLabel = nullptr,
                      const char *outDir    = nullptr)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetPadRightMargin(0.13);

    // ── Paths / labels ────────────────────────────────────────────────────────
    const char *FINE_CSV  =
        "/Users/haggerty/sphenix/data/data02/sphenix/MagnetMapping/cernfinal"
        "/fieldMapFineFullField.csv";
    const char *ROUGH_CSV =
        "/Users/haggerty/sphenix/data/data02/sphenix/MagnetMapping/cernfinal"
        "/fieldMapRoughFullField.csv";
    const char *CALC_ROOT = calcRoot ? calcRoot :
        "fieldmap-used-in-analysis/"
        "8e4d6c3b1660540a658da3a275af2bde_sphenix3dtrackingmapxyz.root";
    const char *LABEL     = calcLabel ? calcLabel : "OPERA";
    const char *OUTDIR    = outDir    ? outDir    : "plots";
    gSystem->mkdir(OUTDIR, kTRUE);

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
    printf("=== Evaluating OPERA map on %d r x %d phi x %d z points ===\n",
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
    double divB_calc [NR][NZ] = {};   // |∇·B|   analysis,  azimuthal mean [T/cm]
    double curlB_calc[NR][NZ] = {};   // |∇×B|   analysis,  azimuthal mean [T/cm]
    double divB_meas [NR][NZ] = {};   // |∇·B|   measured,  cylindrical    [T/cm]
    double curlB_meas[NR][NZ] = {};   // |(∇×B)φ| measured, cylindrical    [T/cm]

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
    // 5b. Maxwell residuals
    //
    // Analysis map: |∇·B| and |∇×B| via central differences (step = DX_C = 2 cm),
    //   evaluated on the comparison cylindrical grid and averaged over phi.
    //
    // Measured map: cylindrical identities on the internal (r,z) grid.
    //   ∇·B  = (1/r)∂(rBr)/∂r + ∂Bz/∂z
    //   (∇×B)φ = ∂Br/∂z − ∂Bz/∂r
    //   Grid steps DR_MM=25 mm, DZ_MM=20 mm → convert to cm for T/cm units.
    // ─────────────────────────────────────────────────────────────────────────
    printf("=== Computing Maxwell residuals ===\n");

    // ── Analysis map ──────────────────────────────────────────────────────────
    {
        const float h     = DX_C;          // 2 cm = 1 grid step
        const float inv2h = 1.f/(2.f*h);
        for (int ir = 0; ir < NR; ++ir) {
            float r = (float)(ir * DR);
            for (int ip = 0; ip < NP; ++ip) {
                double phi = ip * (TMath::TwoPi()/NP);
                float  xq  = r*(float)std::cos(phi);
                float  yq  = r*(float)std::sin(phi);
                for (int iz = 0; iz < NZ; ++iz) {
                    float zq = (float)(Z0 + iz*DZ);
                    float Bx_xp,By_xp,Bz_xp, Bx_xm,By_xm,Bz_xm;
                    float Bx_yp,By_yp,Bz_yp, Bx_ym,By_ym,Bz_ym;
                    float Bx_zp,By_zp,Bz_zp, Bx_zm,By_zm,Bz_zm;
                    CalcInterp(xq+h,yq,  zq,  Bx_xp,By_xp,Bz_xp);
                    CalcInterp(xq-h,yq,  zq,  Bx_xm,By_xm,Bz_xm);
                    CalcInterp(xq,  yq+h,zq,  Bx_yp,By_yp,Bz_yp);
                    CalcInterp(xq,  yq-h,zq,  Bx_ym,By_ym,Bz_ym);
                    CalcInterp(xq,  yq,  zq+h,Bx_zp,By_zp,Bz_zp);
                    CalcInterp(xq,  yq,  zq-h,Bx_zm,By_zm,Bz_zm);
                    double divB = (double(Bx_xp-Bx_xm)+double(By_yp-By_ym)
                                  +double(Bz_zp-Bz_zm))*inv2h;
                    double cBx  = (double(Bz_yp-Bz_ym)-double(By_zp-By_zm))*inv2h;
                    double cBy  = (double(Bx_zp-Bx_zm)-double(Bz_xp-Bz_xm))*inv2h;
                    double cBz  = (double(By_xp-By_xm)-double(Bx_yp-Bx_ym))*inv2h;
                    divB_calc [ir][iz] += std::abs(divB);
                    curlB_calc[ir][iz] += std::sqrt(cBx*cBx+cBy*cBy+cBz*cBz);
                }
            }
        }
        for (int ir = 0; ir < NR; ++ir)
            for (int iz = 0; iz < NZ; ++iz) {
                divB_calc [ir][iz] /= NP;
                curlB_calc[ir][iz] /= NP;
            }
    }

    // ── Measured map (cylindrical, internal grid) ─────────────────────────────
    {
        const double DR_cm  = DR_MM * 0.1;    // 25 mm → 2.5 cm
        const double DZ_cm  = DZ_MM * 0.1;    // 20 mm → 2.0 cm
        const double inv2DZ = 1.0/(2.0*DZ_cm);
        for (int ir = 0; ir < NR; ++ir) {
            int    ir_m = ir;              // 1:1 mapping (same 25-mm step)
            double r_cm = ir_m * DR_cm;
            for (int iz = 0; iz < NZ; ++iz) {
                int iz_m = IZ0 + iz;
                // iz_m in [80,190]: always interior, guards below are safety only
                if (iz_m < 1 || iz_m >= NZ_M-1) continue;

                double dBz_dz = (bz_m[ir_m][iz_m+1]-bz_m[ir_m][iz_m-1])*inv2DZ;

                double dRBr_dr_over_r;
                if (ir_m == 0) {
                    // On axis Br=0; L'Hôpital: (1/r)∂(rBr)/∂r → 2·∂Br/∂r
                    dRBr_dr_over_r = 2.0*br_m[1][iz_m]/DR_cm;
                } else {
                    double rp = (ir_m+1)*DR_cm, rm_r = (ir_m-1)*DR_cm;
                    dRBr_dr_over_r = (rp*br_m[ir_m+1][iz_m]-rm_r*br_m[ir_m-1][iz_m])
                                     /(r_cm*2.0*DR_cm);
                }
                divB_meas[ir][iz] = std::abs(dRBr_dr_over_r+dBz_dz);

                double dBr_dz = (br_m[ir_m][iz_m+1]-br_m[ir_m][iz_m-1])*inv2DZ;
                double dBz_dr = (ir_m > 0 && ir_m < NR_M-1)
                    ? (bz_m[ir_m+1][iz_m]-bz_m[ir_m-1][iz_m])/(2.0*DR_cm) : 0.;
                curlB_meas[ir][iz] = std::abs(dBr_dz-dBz_dr);
            }
        }
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
    TH2F *hBz_calc = MH2("hBz_calc", Form("%s Bz (#phi-avg);z (cm);r (cm)", LABEL));
    TH2F *hdBz     = MH2("hdBz",     Form("#DeltaBz (%s #minus Measured) [mT];z (cm);r (cm)", LABEL));
    TH2F *hBr_meas = MH2("hBr_meas", "Measured Br;z (cm);r (cm)");
    TH2F *hBr_calc = MH2("hBr_calc", Form("%s Br (#phi-avg);z (cm);r (cm)", LABEL));
    TH2F *hdBr     = MH2("hdBr",     Form("#DeltaBr (%s #minus Measured) [mT];z (cm);r (cm)", LABEL));
    TH2F *hBp_m0   = MH2("hBp_m0",  "B#phi  #phi-mean [mT];z (cm);r (cm)");
    TH2F *hBp_A1   = MH2("hBp_A1",  "B#phi  m=1 amplitude [mT];z (cm);r (cm)");
    TH2F *hBz_A1h  = MH2("hBz_A1",  "Bz  m=1 amplitude [mT];z (cm);r (cm)");
    TH2F *hBr_A1h  = MH2("hBr_A1",  "Br  m=1 amplitude [mT];z (cm);r (cm)");
    TH2F *hBp_ph   = MH2("hBp_ph",  "B#phi  m=1 phase [deg];z (cm);r (cm)");
    TH2F *hBz_ph_h    = MH2("hBz_ph",      "Bz  m=1 phase [deg];z (cm);r (cm)");
    TH2F *hDivB_calc  = MH2("hDivB_calc",  Form("%s |#nablaB| [mT/cm];z (cm);r (cm)", LABEL));
    TH2F *hCurlB_calc = MH2("hCurlB_calc", Form("%s |#nabla#timesB| [mT/cm];z (cm);r (cm)", LABEL));
    TH2F *hDivB_meas  = MH2("hDivB_meas",  "Measured |#nablaB| [mT/cm];z (cm);r (cm)");
    TH2F *hCurlB_meas = MH2("hCurlB_meas", "Measured |(#nabla#timesB)_{#phi}| [mT/cm];z (cm);r (cm)");

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
            hBz_ph_h   ->SetBinContent(bx, by, Bz_ph[ir][iz]);
            hDivB_calc ->SetBinContent(bx, by, divB_calc [ir][iz]*1e3);
            hCurlB_calc->SetBinContent(bx, by, curlB_calc[ir][iz]*1e3);
            hDivB_meas ->SetBinContent(bx, by, divB_meas [ir][iz]*1e3);
            hCurlB_meas->SetBinContent(bx, by, curlB_meas[ir][iz]*1e3);
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
        leg->AddEntry(gC, Form("%s (#phi-avg)", LABEL), "l");
        leg->Draw();

        c->cd(2);
        gD->SetTitle(Form("On-axis #DeltaBz = %s #minus Measured;z (cm);#DeltaBz (mT)", LABEL));
        gD->Draw("AL");
        TLine *zl = new TLine(Z0, 0., ZMAX, 0.);
        zl->SetLineStyle(2); zl->SetLineColor(kGray+1); zl->Draw();

        SaveClose(c, Form("%s/01_onaxis_Bz.pdf", OUTDIR));
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
        SaveClose(c, Form("%s/02_Bz_2d_maps.pdf", OUTDIR));
    }

    // ── 03. ΔBz m=0 ───────────────────────────────────────────────────────────
    {
        UseBWR();
        SymZ(hdBz);
        TCanvas *c = new TCanvas("c03", "", 950, 500);
        hdBz->Draw("COLZ");
        SaveClose(c, Form("%s/03_dBz_m0.pdf", OUTDIR));
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
        SaveClose(c, Form("%s/04_Br_comparison.pdf", OUTDIR));
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
        SaveClose(c, Form("%s/05_m1_amplitudes.pdf", OUTDIR));
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

        SaveClose(c, Form("%s/06_m1_amplitude_and_phase.pdf", OUTDIR));
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
        SaveClose(c, Form("%s/07_Bphi_overview.pdf", OUTDIR));
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
                leg->AddEntry(gCr, Form("%s (#phi-avg)", LABEL), "l");
                leg->Draw();
            }
        }
        SaveClose(c, Form("%s/08_Bz_radial_profiles.pdf", OUTDIR));
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
            // Expand y-axis to include measVal if it falls outside the auto-scaled range
            TH1F *axh = g->GetHistogram();
            double ylo = axh->GetMinimum(), yhi = axh->GetMaximum();
            double ymargin = 0.15 * (yhi - ylo);
            if (measVal < ylo) axh->SetMinimum(measVal - ymargin);
            if (measVal > yhi) axh->SetMaximum(measVal + ymargin);
            g->Draw("AL");  // redraw with updated axis
            TLine *ml = new TLine(0., measVal, 350., measVal);
            ml->SetLineColor(kBlue); ml->SetLineWidth(2); ml->SetLineStyle(2);
            ml->Draw();
            TLegend *leg = new TLegend(0.58, 0.77, 0.98, 0.92);
            leg->AddEntry(g,  LABEL,                  "l");
            leg->AddEntry(ml, "Measured (#phi-avg)", "l");
            leg->Draw();
        };
        DrawPhi(1, gBzP, bz_m[ir_s][IZ0+iz_s], "Bz (T)");
        DrawPhi(2, gBrP, br_m[ir_s][IZ0+iz_s], "Br (T)");
        DrawPhi(3, gBpP, 0.,                    "B#phi (T)");
        SaveClose(c, Form("%s/09_phi_dependence_r30_z0.pdf", OUTDIR));
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
        SaveClose(c, Form("%s/10_dBz_phi0_vs_phi180.pdf", OUTDIR));
        delete hdBz0; delete hdBz180;
    }

    // ── 11. Maxwell residuals ─────────────────────────────────────────────────
    // Top row: OPERA map |∇·B| and |∇×B| (both should be ~0 for OPERA FEM solution)
    // Bottom row: measured map |∇·B| (zero by construction) and |(∇×B)_φ| (not enforced)
    {
        gStyle->SetPalette(kViridis);
        TCanvas *c = new TCanvas("c11", "", 1400, 900);
        c->Divide(2, 2);
        c->cd(1); hDivB_calc ->Draw("COLZ");
        c->cd(2); hCurlB_calc->Draw("COLZ");
        c->cd(3); hDivB_meas ->Draw("COLZ");
        c->cd(4); hCurlB_meas->Draw("COLZ");
        SaveClose(c, Form("%s/11_maxwell_residuals.pdf", OUTDIR));
    }

    // ─────────────────────────────────────────────────────────────────────────
    // 8. Save all histograms to a ROOT file
    // ─────────────────────────────────────────────────────────────────────────
    TFile *fOut = TFile::Open(Form("%s/comparison_histograms.root", OUTDIR), "RECREATE");
    hBz_meas->Write(); hBz_calc->Write(); hdBz->Write();
    hBr_meas->Write(); hBr_calc->Write(); hdBr->Write();
    hBp_m0->Write();   hBp_A1->Write();   hBp_ph->Write();
    hBz_A1h->Write();  hBz_ph_h->Write();
    hBr_A1h->Write();
    hDivB_calc->Write(); hCurlB_calc->Write();
    hDivB_meas->Write(); hCurlB_meas->Write();
    fOut->Close();
    printf("  -> %s/", OUTDIR); printf("comparison_histograms.root\n");

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

    // Maxwell residual statistics (r > 0 to skip trivially-zero axis)
    {
        double maxDivC=0., maxCurlC=0., maxDivM=0., maxCurlM=0.;
        double rmssDivC=0., rmssCurlC=0., rmssDivM=0., rmssCurlM=0.;
        int n = 0;
        for (int ir = 1; ir < NR; ++ir)
            for (int iz = 0; iz < NZ; ++iz) {
                maxDivC  = std::max(maxDivC,  divB_calc [ir][iz]);
                maxCurlC = std::max(maxCurlC, curlB_calc[ir][iz]);
                maxDivM  = std::max(maxDivM,  divB_meas [ir][iz]);
                maxCurlM = std::max(maxCurlM, curlB_meas[ir][iz]);
                rmssDivC  += divB_calc [ir][iz]*divB_calc [ir][iz];
                rmssCurlC += curlB_calc[ir][iz]*curlB_calc[ir][iz];
                rmssDivM  += divB_meas [ir][iz]*divB_meas [ir][iz];
                rmssCurlM += curlB_meas[ir][iz]*curlB_meas[ir][iz];
                ++n;
            }
        double sc = 1e3/std::sqrt((double)n);   // T/cm → mT/cm, /sqrt(n)
        printf("\n  Maxwell residuals (r > 0, tracking volume):\n");
        printf("  %-9s |div B|:        max = %7.3f mT/cm,  RMS = %7.4f mT/cm\n",
               LABEL, maxDivC*1e3,  std::sqrt(rmssDivC )*sc);
        printf("  %-9s |curl B|:       max = %7.3f mT/cm,  RMS = %7.4f mT/cm\n",
               LABEL, maxCurlC*1e3, std::sqrt(rmssCurlC)*sc);
        printf("  Measured  |div B|:        max = %7.3f mT/cm,  RMS = %7.4f mT/cm\n",
               maxDivM*1e3,  std::sqrt(rmssDivM )*sc);
        printf("  Measured  |(curl B)_phi|: max = %7.3f mT/cm,  RMS = %7.4f mT/cm\n",
               maxCurlM*1e3, std::sqrt(rmssCurlM)*sc);
        printf("  (Analysis step = %.0f cm; Measured step DR=%.1f cm, DZ=%.1f cm)\n",
               (double)DX_C, DR_MM*0.1, DZ_MM*0.1);
    }

    // Clean up
    delete[] gCalcBx; delete[] gCalcBy; delete[] gCalcBz;
    gCalcBx = gCalcBy = gCalcBz = nullptr;
    printf("\nDone.\n");
}
