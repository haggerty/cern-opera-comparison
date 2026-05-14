// checkTilt.C
//
// Azimuthal Fourier analysis of the raw CERN measurement CSVs.
// Searches for the m=1 dipole component of Bz(r,phi,z).
//
// A tilt of angle alpha toward azimuth phi0 produces:
//   dBz(r,phi,z) ~ -(dBz/dr) * z * alpha * cos(phi - phi0)
// The amplitude is proportional to |z| and reverses sign at z=0,
// distinguishing a tilt from a simple axis translation.
//
// For a 2.39 mrad survey tilt toward +y (phi0 ~ +90 deg):
//   expected m=1 phase ~ +90 deg for z > 0 where dBz/dr < 0
//
// CSV format (surveyor frame): x_s,y_s,z_s,|B|,Bx_s,By_s,Bz_s  [mm, T]
// sPHENIX transform: xp=xs, yp=zs, zp=-ys; Bz_phx = -By_s
//
// Run: root -l -b -q 'checkTilt.C+'

#include <cstdio>
#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <algorithm>

#include "TCanvas.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"

// Grid matching sPHENIXFieldMap
static const int    kNR   = 37,    kNZ   = 241;
static const double kRMin =   0.,  kRMax = 900., kdR =  25.;  // mm
static const double kZMin = -2700.,kZMax = 2100., kdZ =  20.;  // mm

struct Meas { double phi, Bz; };
using Key = std::pair<int,int>;

// Least-squares fit: Bz(phi) = A0 + A1*cos(phi) + B1*sin(phi)
static bool FitFourier(const std::vector<Meas>& v,
                       double& A0, double& A1, double& B1)
{
    const int N = (int)v.size();
    if (N < 3) return false;
    double S1=0,Sc=0,Ss=0,Sc2=0,Ss2=0,Scs=0,SY=0,SYc=0,SYs=0;
    for (const auto& m : v) {
        double c=cos(m.phi), s=sin(m.phi);
        S1+=1; Sc+=c; Ss+=s; Sc2+=c*c; Ss2+=s*s; Scs+=c*s;
        SY+=m.Bz; SYc+=m.Bz*c; SYs+=m.Bz*s;
    }
    double a[3][4]={{S1,Sc,Ss,SY},{Sc,Sc2,Scs,SYc},{Ss,Scs,Ss2,SYs}};
    for (int col=0; col<3; ++col) {
        int piv=col;
        for (int row=col+1; row<3; ++row)
            if (fabs(a[row][col]) > fabs(a[piv][col])) piv=row;
        if (fabs(a[piv][col]) < 1e-15) return false;
        if (piv!=col)
            for (int j=0; j<4; ++j) std::swap(a[col][j], a[piv][j]);
        double fac=a[col][col];
        for (int j=col; j<4; ++j) a[col][j]/=fac;
        for (int row=0; row<3; ++row) {
            if (row==col) continue;
            double f=a[row][col];
            for (int j=col; j<4; ++j) a[row][j]-=f*a[col][j];
        }
    }
    A0=a[0][3]; A1=a[1][3]; B1=a[2][3];
    return true;
}

static void ReadCSV(const char* fname, bool isFine,
                    std::map<Key, std::vector<Meas>>& bins,
                    std::set<Key>& fineBins)
{
    FILE* f = fopen(fname, "r");
    if (!f) { fprintf(stderr,"Cannot open %s\n",fname); return; }
    printf("Reading %s ...\n", fname);
    long nread=0, nused=0;
    char line[512];
    while (fgets(line, sizeof(line), f)) {
        if (line[0]=='#' || line[0]=='\n') continue;
        double xs,ys,zs,Bmag,Bxs,Bys,Bzs;
        if (sscanf(line,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   &xs,&ys,&zs,&Bmag,&Bxs,&Bys,&Bzs)!=7) continue;
        ++nread;
        double xp=xs, yp=zs, zp=-ys;
        double r  =sqrt(xp*xp+yp*yp);
        double phi=atan2(yp,xp);
        double Bz =-Bys;
        if (r <kRMin || r >kRMax+kdR) continue;
        if (zp<kZMin || zp>kZMax+kdZ) continue;
        int ir=(int)round((r -kRMin)/kdR);
        int iz=(int)round((zp-kZMin)/kdZ);
        if (ir<0) ir=0; if (ir>=kNR) ir=kNR-1;
        if (iz<0) iz=0; if (iz>=kNZ) iz=kNZ-1;
        Key key={ir,iz};
        if (!isFine && fineBins.count(key)) continue;
        bins[key].push_back({phi,Bz});
        if (isFine) fineBins.insert(key);
        ++nused;
    }
    fclose(f);
    printf("  %ld rows read, %ld used\n", nread, nused);
}

static void SaveClose(TCanvas* c, const char* path) { c->SaveAs(path); delete c; }

// ─────────────────────────────────────────────────────────────────────────────
void checkTilt(
    const char* fineCSV  = "/Users/haggerty/sphenix/data/data02/sphenix/MagnetMapping/cernfinal/fieldMapFineFullField.csv",
    const char* roughCSV = "/Users/haggerty/sphenix/data/data02/sphenix/MagnetMapping/cernfinal/fieldMapRoughFullField.csv",
    const char* outDir   = "plots_tilt")
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    // ── Load data ─────────────────────────────────────────────────────────────
    std::map<Key, std::vector<Meas>> bins;
    std::set<Key> fineBins;
    ReadCSV(fineCSV,  true,  bins, fineBins);
    ReadCSV(roughCSV, false, bins, fineBins);
    printf("Total: %zu (ir,iz) bins\n\n", bins.size());

    // Print phi-count distribution
    {
        int cnt[32]={};
        for (auto& kv : bins) {
            int n=std::min((int)kv.second.size(),31);
            cnt[n]++;
        }
        printf("phi values per bin:\n");
        for (int i=0;i<32;++i) if (cnt[i]) printf("  N=%2d: %5d bins\n",i,cnt[i]);
        printf("\n");
    }

    // ── Histograms ─────────────────────────────────────────────────────────────
    // Each bin centered on a grid node; z in cm, r in mm
    const double drH=kdR/2., dzH=kdZ/2.;
    const double zlo=(kZMin-dzH)/10., zhi=(kZMax+dzH)/10.;
    const double rlo=kRMin-drH,       rhi=kRMax+drH;

    TH2D* hNphi  =new TH2D("hNphi",  "N(#phi) per (r,z) bin;z (cm);r (mm)",
                            kNZ,zlo,zhi, kNR,rlo,rhi);
    TH2D* hA1frac=new TH2D("hA1frac","m=1/m=0 amplitude;z (cm);r (mm)",
                            kNZ,zlo,zhi, kNR,rlo,rhi);
    TH2D* hA1mT  =new TH2D("hA1mT",  "m=1 amplitude (mT);z (cm);r (mm)",
                            kNZ,zlo,zhi, kNR,rlo,rhi);
    TH2D* hPhase =new TH2D("hPhase", "m=1 phase (deg);z (cm);r (mm)",
                            kNZ,zlo,zhi, kNR,rlo,rhi);

    // Profiles at selected r nodes (r = 50, 100, 200, 300, 450, 600, 750, 800 mm)
    const int selIR[]={2,4,8,12,18,24,30,32};
    const int nSel=8;
    int colors[]={kBlue+1,kRed+1,kGreen+2,kMagenta+1,kOrange+2,kCyan+2,kViolet+1,kGray+2};
    TProfile *pPhZ[nSel], *pA1Z[nSel];
    for (int k=0;k<nSel;++k) {
        int ir=selIR[k];
        double rmm=kRMin+ir*kdR;
        pPhZ[k]=new TProfile(Form("pPhZ_%d",ir),
            Form("r=%.0f mm;z (cm);m=1 phase (deg)",rmm),
            kNZ,zlo,zhi,-200.,200.);
        pA1Z[k]=new TProfile(Form("pA1Z_%d",ir),
            Form("r=%.0f mm;z (cm);A_{1}/A_{0}",rmm),
            kNZ,zlo,zhi,-2.,2.);
        for (auto* p:{pPhZ[k],pA1Z[k]}) {
            p->SetLineColor(colors[k]);
            p->SetMarkerColor(colors[k]);
            p->SetMarkerStyle(20);
            p->SetMarkerSize(0.8);
        }
    }

    // ── Fourier loop ───────────────────────────────────────────────────────────
    int nFit=0, nSkip=0;
    double phSum=0.; int phN=0;
    for (auto& kv : bins) {
        int ir=kv.first.first, iz=kv.first.second;
        auto& v=kv.second;
        double z_cm=(kZMin+iz*kdZ)/10.;
        double r_mm=kRMin+ir*kdR;
        hNphi->Fill(z_cm, r_mm, (double)v.size());
        double A0,A1,B1;
        if (!FitFourier(v,A0,A1,B1)){++nSkip;continue;}
        ++nFit;
        double amp=sqrt(A1*A1+B1*B1);
        double ph =atan2(B1,A1)*(180./TMath::Pi());
        double frac=(fabs(A0)>1e-6) ? amp/fabs(A0) : 0.;
        hA1frac->Fill(z_cm,r_mm,frac);
        hA1mT  ->Fill(z_cm,r_mm,amp*1e3);
        hPhase ->Fill(z_cm,r_mm,ph);
        // Accumulate mean phase in central tracking region
        if (fabs(z_cm)<150. && r_mm>0. && r_mm<600. && A0>0.5)
            {phSum+=ph; ++phN;}
        for (int k=0;k<nSel;++k) {
            if (ir==selIR[k]) {
                pPhZ[k]->Fill(z_cm,ph);
                if (A0>0.5) pA1Z[k]->Fill(z_cm,frac);
            }
        }
    }
    printf("Fourier: %d fits, %d skipped (<3 phi values)\n", nFit, nSkip);
    if (phN>0)
        printf("Mean m=1 phase (|z|<150cm, r<600mm, A0>0.5T): %.1f deg  [N=%d]\n",
               phSum/phN, phN);
    printf("  Survey direction: phi0 ~ +71 deg\n");
    printf("  OPERA m=1 phase: -86 deg\n\n");

    // ── Plots ──────────────────────────────────────────────────────────────────
    gSystem->mkdir(outDir, kTRUE);

    // A: N(phi) per bin — shows measurement coverage
    {
        auto* c=new TCanvas("cNphi","",900,450);
        hNphi->Draw("colz");
        SaveClose(c, Form("%s/tilt_A_nphi_per_bin.pdf",outDir));
    }
    // B: m=1/m=0 fractional amplitude map
    {
        auto* c=new TCanvas("cA1f","",900,450);
        hA1frac->Draw("colz");
        SaveClose(c, Form("%s/tilt_B_m1_frac_map.pdf",outDir));
    }
    // C: m=1 absolute amplitude map
    {
        auto* c=new TCanvas("cA1mT","",900,450);
        hA1mT->Draw("colz");
        SaveClose(c, Form("%s/tilt_C_m1_abs_map.pdf",outDir));
    }
    // D: m=1 phase map
    {
        auto* c=new TCanvas("cPhMap","",900,450);
        hPhase->Draw("colz");
        SaveClose(c, Form("%s/tilt_D_phase_map.pdf",outDir));
    }
    // E: phase vs z at selected r values
    {
        auto* c=new TCanvas("cPhZ","",900,600);
        c->SetGrid();
        gPad->SetTopMargin(0.24);
        gPad->SetBottomMargin(0.12);
        bool first=true;
        for (int k=0;k<nSel;++k) {
            if (pPhZ[k]->GetEntries()<1) continue;
            pPhZ[k]->SetTitle(";z (cm);m=1 phase (deg)");
            pPhZ[k]->GetYaxis()->SetRangeUser(-200.,200.);
            pPhZ[k]->Draw(first ? "E" : "E SAME");
            first=false;
        }
        // Reference lines
        auto* lSurv=new TLine(-270.,71.,210.,71.);
        auto* lOPER=new TLine(-270.,-86.,210.,-86.);
        lSurv->SetLineStyle(2); lSurv->SetLineColor(kCyan+2);   lSurv->Draw();
        lOPER->SetLineStyle(2); lOPER->SetLineColor(kOrange+1); lOPER->Draw();
        TLatex lt; lt.SetTextSize(0.033);
        lt.SetTextColor(kCyan+2);  lt.DrawLatex(-260., 80.,"Survey: #phi_{0}=+71#circ");
        lt.SetTextColor(kOrange+1);lt.DrawLatex(-260.,-98.,"OPERA:  #phi_{0}=#minus86#circ");
        // Title and legend in top margin
        TLatex tt; tt.SetNDC(); tt.SetTextAlign(22);
        tt.SetTextSize(0.042); tt.SetTextFont(42);
        tt.DrawLatex(0.5, 0.945, "Measured map: m=1 Bz phase vs z");
        TLegend* leg=new TLegend(0.08,0.76,0.92,0.92);
        leg->SetNColumns(4);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.038);
        for (int k=0;k<nSel;++k)
            if (pPhZ[k]->GetEntries()>0)
                leg->AddEntry(pPhZ[k],Form("r = %.0f mm",kRMin+selIR[k]*kdR),"lp");
        leg->Draw();
        SaveClose(c, Form("%s/tilt_E_phase_vs_z.pdf",outDir));
    }
    // F: amplitude vs z at selected r values
    {
        auto* c=new TCanvas("cA1Z","",900,550);
        c->SetGrid();
        TLegend* leg=new TLegend(0.65,0.65,0.88,0.88);
        bool first=true;
        for (int k=0;k<nSel;++k) {
            if (pA1Z[k]->GetEntries()<1) continue;
            pA1Z[k]->SetTitle("Measured map: m=1/m=0 amplitude vs z");
            if (first) pA1Z[k]->SetMaximum(0.015);
            pA1Z[k]->Draw(first ? "E" : "E SAME");
            leg->AddEntry(pA1Z[k],Form("r = %.0f mm",kRMin+selIR[k]*kdR),"lp");
            first=false;
        }
        leg->Draw();
        SaveClose(c, Form("%s/tilt_F_amplitude_vs_z.pdf",outDir));
    }
    printf("Plots saved to %s/tilt_*.pdf\n", outDir);
}
