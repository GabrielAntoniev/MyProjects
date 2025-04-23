#include <iostream>
#include <graphics.h>
#include <cmath>
#include <fstream>
#include <windowsx.h>
using namespace std;

struct plan{double a,b,c,d;};
struct pereche{int x,y;};
struct punct{double x,y,z;};
struct muchie{punct p1, p2; muchie* urm;};muchie *prim=NULL, *ultim=NULL;
struct puncteCG{punct treiD; pereche doiD; muchie *inceput, *sfarsit; puncteCG *urm;};puncteCG *primCG=NULL, *ultimCG=NULL;
///struct puncte3D{punct treiD; pereche doiD; muchie *poz; puncte3D *urm;};puncte3D *prim3D=NULL, *ultim3D=NULL;

int culoare_linie, fundal = COLOR(195, 217, 247);
int page;
punct eye, origine, centru_window, c1, c2, c3, c4, d1, d2, d3, d4;
plan window;
double dist_eye_centru_window;
double pi=3.14159265358979323846264, doipi=pi+pi;
bool tastare = false, alt_window_open = false, alt_window_open2 = false, anulare = false, intro = true;
bool romana_selected = true, light_selected = true;
int color = 0, fundal_intro = COLOR(195, 217, 247), culoare_butoane;

int N=1050, M=795;
int Npe2=N/2, Mpe2=M/2;
int NRMUCHII;
int pagForme, pagPrincipala, pagCrt, pagSave, pagOpen, pagIntro;
double modul(double x){if(x>0)return x;return -x;}


plan ecuatie_plan(punct p1, punct p2, punct p3);
void initializare();
punct simetric(punct p, punct s1, punct s2);
punct proiectie(punct p, punct s1, punct s2);
double dist_punct_punct(punct p1, punct p2);
double dist_punct_dreapta(punct p, punct p1, punct p2);
double dist_punct_plan(punct p, plan q);
punct intersectie_plan(punct p);
pereche choose_pixel(punct p);
///bool se_vede(punct p);
void proiectare(bool start);
void ajustare_puncte_principale_la_zoom(double d, bool b);
punct mijloc_segment(punct p1, punct p2);
void ajustare_puncte_principale_la_rotatie(pereche p1, pereche p2);
punct rotatie(punct p, double unghi, punct u);
punct normalizare(punct r);
puncteCG* gasitCG(pereche s);
void translare_corp(pereche p1, pereche p2, puncteCG* adresa);
void rotire_corp(pereche p1, pereche p2, puncteCG* &adr);
void translare_punct(pereche p1, pereche p2, muchie* adr, punct &pun);
void proiectareCG();
void stergere();
void rutina_test(bool start_to_project_from_beggining);
void gui_forme_menu();
void gui_forme_shapes();
void gui_save();
void adaugare_fisier_obj(punct c, int L, char numeFisier[]);
void adaugare_fisier_txt(punct c, int L, char numeFisier[]);
/////////////////////////////////////////////
punct rotatie_y(punct p);
/////////////////////////////////////////////////


plan ecuatie_plan(punct p1, punct p2, punct p3)
{
    ///det ecuatie plan stiind 3 puncte 3D
    plan rez;
    double yp21=p2.y-p1.y, zp31=p3.z-p1.z, yp31=p3.y-p1.y, zp21=p2.z-p1.z, xp21=p2.x-p1.x, xp31=p3.x-p1.x;

    rez.a=(yp21)*(zp31)-(yp31)*(zp21);
    rez.b=(xp31)*(zp21)-(xp21)*(zp31);
    rez.c=(xp21)*(yp31)-(yp21)*(xp31);
    rez.d=-rez.a*p1.x-rez.b*p1.y-rez.c*p1.z;
    return rez;
}

void initializare()
{
    eye.x=0.0;eye.y=0.0;eye.z=2000.0;origine.x=0.0;origine.y=0.0;origine.z=0.0;
    centru_window.x=0.0;centru_window.y=0.0;centru_window.z=1700.0;

    dist_eye_centru_window=eye.z-centru_window.z;

    c1.x=-N/2; c1.y=M/2; c1.z=centru_window.z;
    c2.x=N/2; c2.y=M/2; c2.z=centru_window.z;
    c3.x=N/2; c3.y=-M/2; c3.z=centru_window.z;
    c4.x=-N/2; c4.y=-M/2; c4.z=centru_window.z;

    d1.x=c1.x; d1.y=0; d1.z=c1.z;
    d2.x=0; d2.y=c1.y; d2.z=c1.z;
    d3.x=c2.x; d3.y=0; d3.z=c1.z;
    d4.x=0; d4.y=c4.y; d4.z=c1.z;

    window=ecuatie_plan(c1,c2,c3);

    prim=new muchie;ultim=prim;
    punct k1=origine; k1.y=100;
    ultim->p1=origine; ultim->p2=k1;

    muchie*t=new muchie;
    k1=origine; k1.z=100;
    ultim->urm=t; ultim=t;
    ultim->p1=origine; ultim->p2=k1;

    t=new muchie;
    k1=origine; k1.x=100;
    ultim->urm=t; ultim=t;
    ultim->p1=origine; ultim->p2=k1; ultim->urm=NULL;
}

void ajustare_puncte_principale_la_zoom(double d, bool b)
{
    ///zoom in : 1
    ///zoom out : 0
    double K_fata_sau_spate=0;
    double dist_eye_origine=dist_punct_punct(eye,origine);

    if(b==true) K_fata_sau_spate=1-d/dist_eye_origine;///folosim fata
    else K_fata_sau_spate=1+d/dist_eye_origine;

        eye.x*=K_fata_sau_spate;eye.y*=K_fata_sau_spate;eye.z*=K_fata_sau_spate;
        punct centru_window_vechi=centru_window;
        double LAMBDA=(1.0)/(1+dist_eye_centru_window/(dist_punct_punct(eye,origine)-dist_eye_centru_window));
        centru_window.x=eye.x*LAMBDA; centru_window.y=eye.y*LAMBDA; centru_window.z=eye.z*LAMBDA;

        c1.x=centru_window.x-centru_window_vechi.x+c1.x;
        c1.y=centru_window.y-centru_window_vechi.y+c1.y;
        c1.z=centru_window.z-centru_window_vechi.z+c1.z;

        c2.x=centru_window.x-centru_window_vechi.x+c2.x;
        c2.y=centru_window.y-centru_window_vechi.y+c2.y;
        c2.z=centru_window.z-centru_window_vechi.z+c2.z;

        c3.x=centru_window.x+centru_window.x-c1.x;
        c3.y=centru_window.y+centru_window.y-c1.y;
        c3.z=centru_window.z+centru_window.z-c1.z;

        c4.x=centru_window.x+centru_window.x-c2.x;
        c4.y=centru_window.y+centru_window.y-c2.y;
        c4.z=centru_window.z+centru_window.z-c2.z;

    d1=mijloc_segment(c1,c4);
    d2=mijloc_segment(c1,c2);
    d3=mijloc_segment(c2,c3);
    d4=mijloc_segment(c3,c4);

    //window=ecuatie_plan(c1,c2,c3);
}

punct normalizare(punct r)
{
    double D=(1.0)/dist_punct_punct(r,origine);
    r.x*=D;r.y*=D;r.z*=D;
    return r;
}

punct rotatie_y(punct p)
{
    double cost, sint;
    cost=cos(pi/90);sint=sin(pi/90);
    punct rez={cost*p.x-sint*p.y, sint*p.x+cost*p.y, p.z};
    return rez;
}

punct rotatie(punct p, double unghi, punct u)
{
   punct rez={0,0,0};
   double cosunghi=cos(unghi), sinunghi=sin(unghi), unuminuscosunghi=1-cosunghi;

   u=normalizare(u);
   double uxuy=u.x*u.y, uyuz=u.y*u.z, uxuz=u.x*u.z;

   rez.x=(cosunghi+unuminuscosunghi*u.x*u.x)*p.x
   +(unuminuscosunghi*uxuy-u.z*sinunghi)*p.y
   +(unuminuscosunghi*uxuz+u.y*sinunghi)*p.z;

   rez.y=(unuminuscosunghi*uxuy+u.z*sinunghi)*p.x
   +(cosunghi+unuminuscosunghi*u.y*u.y)*p.y
   +(unuminuscosunghi*uyuz-u.x*sinunghi)*p.z;

   rez.z=(unuminuscosunghi*uxuz-u.y*sinunghi)*p.x
   +(unuminuscosunghi*uyuz+u.x*sinunghi)*p.y
   +(cosunghi+unuminuscosunghi*u.z*u.z)*p.z;

   return rez;
}

void ajustare_puncte_principale_la_rotatie(pereche p1, pereche p2)
{
    ///p1 si p2 reprezinta vectorul contruit cand dau Rclick si trag mouse ul pe ecran
    double L, K, distx, disty, U=0.000005;
    punct MIJ;

    if((p2.x>p1.x&&p1.y>=p2.y)||(p1.x>p2.x&&p2.y>=p1.y))///caz 1 si 4
    {
        if(p1.x>p2.x&&p2.y>=p1.y)U=-0.000005;

        distx=modul(p1.x-p2.x);disty=modul(p1.y-p2.y);
        L=disty/(Mpe2-disty);K=(1.0)/(1+L);
        punct k1={(centru_window.x+L*d2.x)*K, (centru_window.y+L*d2.y)*K, (centru_window.z+L*d2.z)*K};

        L=distx/(Npe2-distx);K=(1.0)/(1+L);
        punct k2={(centru_window.x+L*d3.x)*K, (centru_window.y+L*d3.y)*K, (centru_window.z+L*d3.z)*K};

        punct mij=mijloc_segment(k1,k2);
        plan P=ecuatie_plan(mij,centru_window,origine);///Planul P este un plan perpendicular pe planul window, cu o unghiulatie data si care se intersecteaza in centru_window

        ///refolosim mij -> devine normala pe plan

        mij.x=P.a; mij.y=P.b; mij.z=P.c; MIJ=mij;
        U*=sqrt((distx*distx+disty*disty)*(eye.x*eye.x+eye.y*eye.y+eye.z*eye.z));

        eye=rotatie(eye,U,MIJ);
        double LAMBDA=(1.0)/(1+dist_eye_centru_window/(dist_punct_punct(eye,origine)-dist_eye_centru_window));
        centru_window.x=eye.x*LAMBDA; centru_window.y=eye.y*LAMBDA; centru_window.z=eye.z*LAMBDA;

        c1=rotatie(c1,U,MIJ);
        c2=rotatie(c2,U,MIJ);
        c3.x=centru_window.x+centru_window.x-c1.x; c3.y=centru_window.y+centru_window.y-c1.y; c3.z=centru_window.z+centru_window.z-c1.z;
        c4.x=centru_window.x+centru_window.x-c2.x; c4.y=centru_window.y+centru_window.y-c2.y; c4.z=centru_window.z+centru_window.z-c2.z;

        d1=mijloc_segment(c1,c4);
        d2=mijloc_segment(c1,c2);
        d3=mijloc_segment(c2,c3);
        d4=mijloc_segment(c3,c4);

        //window=ecuatie_plan(c1,c2,c3);
    }

    else if((p1.x>=p2.x&&p1.y>p2.y)||(p2.x>=p1.x&&p2.y>p1.y))///caz 3 si 2
    {
        if(p2.x>=p1.x&&p2.y>p1.y)U=-0.000005;

        distx=modul(p1.x-p2.x);disty=modul(p1.y-p2.y);
        L=disty/(Mpe2-disty);K=(1.0)/(1+L);
        punct k1={(centru_window.x+L*d2.x)*K, (centru_window.y+L*d2.y)*K, (centru_window.z+L*d2.z)*K};

        L=distx/(Npe2-distx);K=(1.0)/(1+L);
        punct k2={(centru_window.x+L*d1.x)*K, (centru_window.y+L*d1.y)*K, (centru_window.z+L*d1.z)*K};

        punct mij=mijloc_segment(k1,k2);
        plan P=ecuatie_plan(mij,centru_window,origine);///Planul P este un plan perpendicular pe planul window, cu o unghiulatie data si care se intersecteaza in centru_window

        ///refolosim mij -> devine normala pe plan

        mij.x=P.a; mij.y=P.b; mij.z=P.c; MIJ=mij;
        U*=sqrt((distx*distx+disty*disty)*(eye.x*eye.x+eye.y*eye.y+eye.z*eye.z));

        eye=rotatie(eye,U,MIJ);
        double LAMBDA=(1.0)/(1+dist_eye_centru_window/(dist_punct_punct(eye,origine)-dist_eye_centru_window));
        centru_window.x=eye.x*LAMBDA; centru_window.y=eye.y*LAMBDA; centru_window.z=eye.z*LAMBDA;

        c1=rotatie(c1,U,MIJ);
        c2=rotatie(c2,U,MIJ);
        c3.x=centru_window.x+centru_window.x-c1.x; c3.y=centru_window.y+centru_window.y-c1.y; c3.z=centru_window.z+centru_window.z-c1.z;
        c4.x=centru_window.x+centru_window.x-c2.x; c4.y=centru_window.y+centru_window.y-c2.y; c4.z=centru_window.z+centru_window.z-c2.z;

        d1=mijloc_segment(c1,c4);
        d2=mijloc_segment(c1,c2);
        d3=mijloc_segment(c2,c3);
        d4=mijloc_segment(c3,c4);

        //window=ecuatie_plan(c1,c2,c3);
    }
}

punct simetric(punct p, punct s1, punct s2)
{
    ///s1 si s2 repr dreapta fata de care se face simetricul
    punct rez={-p.x,-p.y,-p.z};
    double dx=s2.x-s1.x, dy=s2.y-s1.y, dz=s2.z-s1.z;
    double t=( dx*(p.x-s1.x)+dy*(p.y-s1.y)+dz*(p.z-s1.z) )/( dx*dx + dy*dy + dz*dz );

    rez.x+= 2*(s1.x+t*dx);
    rez.y+= 2*(s1.y+t*dy);
    rez.z+= 2*(s1.z+t*dz);
    return rez;
}

punct proiectie(punct p, punct s1, punct s2)
{
    ///proiectia lui p se face pe dreapta determinata de punctele s1 si s2
    double dx=s2.x-s1.x, dy=s2.y-s1.y, dz=s2.z-s1.z;
    double t=( dx*(p.x-s1.x)+dy*(p.y-s1.y)+dz*(p.z-s1.z) )/( dx*dx + dy*dy + dz*dz );
    punct rez={s1.x+t*dx, s1.y+t*dy, s1.z+t*dz};
    return rez;
}

punct mijloc_segment(punct p1, punct p2)
{
    punct rez={0.5*(p1.x+p2.x), 0.5*(p1.y+p2.y), 0.5*(p1.z+p2.z)};
    return rez;
}

double dist_punct_punct(punct p1, punct p2)
{
    double r1=p1.x-p2.x, r2=p1.y-p2.y, r3=p1.z-p2.z;
    return sqrt(r1*r1+r2*r2+r3*r3);
}

double dist_punct_dreapta(punct p, punct p1, punct p2)
{
    /// dreapta formata de p1 si p2
    ///return dist_punct_punct(p,proiectie(p,p1,p2));
    double r1=p1.x-p2.x, r2=p1.y-p2.y, r3=p1.z-p2.z;

    double x1=p1.x-p.x, x2=p2.x-p.x, y1=p1.y-p.y, y2=p2.y-p.y, z1=p1.z-p.z, z2=p2.z-p.z;
    double c1=y1*z2-z1*y2, c2=z1*x2-x1*z2, c3=x1*y2-y1*x2;
    return sqrt((c1*c1+c2*c2+c3*c3)/(r1*r1+r2*r2+r3*r3));
}

double dist_punct_plan(punct p, plan q)
{
    ///va fi "signed" distance, adica valoarea returnata va fi pozitiva daca punctul p sta pe partea unde normala pe plan este pozitiva
                             ///si distanta va fi negativa daca punctul sta pe partea de plan unde normala planului este negativa


}
/*
pereche choose_pixel_super(punct p)
{
    double r1,r2,r3,r4;
    double dist_c1_patrat,dist_c2_patrat,dist_c3_patrat,dist_c4_patrat;

    r1=p.x-c1.x, r2=p.y-c1.y, r3=p.z-c1.z;
    dist_c1_patrat=r1*r1+r2*r2+r3*r3;

    r1=p.x-c2.x, r2=p.y-c2.y, r3=p.z-c2.z;
    dist_c2_patrat=r1*r1+r2*r2+r3*r3;

    r1=p.x-c3.x, r2=p.y-c3.y, r3=p.z-c3.z;
    dist_c3_patrat=r1*r1+r2*r2+r3*r3;

    r1=p.x-c4.x, r2=p.y-c4.y, r3=p.z-c4.z;
    dist_c4_patrat=r1*r1+r2*r2+r3*r3;

    double dist_min=min(min(min(dist_c1_patrat,dist_c2_patrat),dist_c3_patrat),dist_c4_patrat);

    double v1,v2=0;
    if(dist_min==dist_c1_patrat)v1=-1,v2=1;
    if(dist_min==dist_c2_patrat)v1=1,v2=1;
    if(dist_min==dist_c3_patrat)v1=1,v2=-1;
    v1=-1,v2=-1;

    plan d1d3=ecuatie_plan(d2,d4,origine);
    plan d2d4=ecuatie_plan(d1,d3,origine);
    double dist_plan_d1d3=

}
*/
punct intersectie_plan(punct p)
{
    ///rezultatul este punctul de intersectie al dreptei eye-p cu planul window
    double ei, ep=dist_punct_punct(eye,p);

    double lambda=dist_eye_centru_window*sqrt(eye.x*eye.x+eye.y*eye.y+eye.z*eye.z)/(eye.x*(eye.x-p.x)+eye.y*(eye.y-p.y)+eye.z*(eye.z-p.z));
    ei=lambda*ep;

    double L=ei/(ep-ei); double K=(1.0)/(1+L);
    punct rez={(eye.x+L*p.x)*K, (eye.y+L*p.y)*K, (eye.z+L*p.z)*K};
    return rez;
}

pereche choose_pixel(punct p)
{
    pereche rez;
    double D1=dist_punct_dreapta(p,d2,d4), D2=dist_punct_dreapta(p,d1,d3);
    double Dp14=dist_punct_dreapta(p,c1,c4), Dp12=dist_punct_dreapta(p,c1,c2);

    rez.x=(int)Dp14; rez.y=(int)Dp12;
    if(D1<=Npe2 && D2<=Mpe2)return rez;

    if(D1>Npe2&&Dp14<=D1)rez.x=-rez.x;
    if(D2>Mpe2&&Dp12<=D2)rez.y=-rez.y;

    return rez;
}

void stergere()
{
    if(primCG!=NULL)
    {
        puncteCG* aux=primCG->urm;
        while(aux!=NULL){
            delete primCG;
            primCG=aux;
            aux=aux->urm;
        }
        delete primCG;
    }

    if(prim!=NULL)
    {
        muchie* aux=prim->urm;
        while(aux!=NULL){
            delete prim;
            prim=aux;
            aux=aux->urm;
        }
        delete prim;
    }

    prim=new muchie;ultim=prim;
    punct k1=origine; k1.y=100;
    ultim->p1=origine; ultim->p2=k1;

    muchie*t=new muchie;
    k1=origine; k1.z=100;
    ultim->urm=t; ultim=t;
    ultim->p1=origine; ultim->p2=k1;

    t=new muchie;
    k1=origine; k1.x=100;
    ultim->urm=t; ultim=t;
    ultim->p1=origine; ultim->p2=k1; ultim->urm=NULL;

    primCG=NULL;
    ultimCG=NULL;
    delay(100);
}

void proiectare(bool start)
{
    ///0- incep de la prim->urm->urm->urm
    ///1-incep normal, cu axele alea rosu, verge si albastru

    muchie*t=prim;
    if(start==0 && t->urm->urm->urm!=NULL)t=t->urm->urm->urm;
    double D=eye.x*eye.x+eye.y*eye.y+eye.z*eye.z;
    do
    {
        if((t->p1).x*eye.x+(t->p1).y*eye.y+(t->p1).z*eye.z<=D && (t->p2).x*eye.x+(t->p2).y*eye.y+(t->p2).z*eye.z<=D){

        pereche rez1=choose_pixel( intersectie_plan(t->p1) );
        pereche rez2=choose_pixel( intersectie_plan(t->p2) );

        if(t==prim)setcolor(RED);//y
        else if(t==prim->urm)setcolor(GREEN);//z
        else if(t==prim->urm->urm)setcolor(BLUE);//x
        //else setcolor(BLUE);
        else setcolor(culoare_linie);
        line(rez1.x, rez1.y, rez2.x, rez2.y);
        }
        t=t->urm;

    }while(t!=NULL);
}

void proiectareCG()
{
    puncteCG*t=primCG;
    setcolor(MAGENTA);
    double D=eye.x*eye.x+eye.y*eye.y+eye.z*eye.z;
    if(primCG!=NULL)
    do
    {
        (t->doiD).x=-9999; (t->doiD).y=-9999;
        double X=0,Y=0,Z=0,nr=0;
        muchie*m=t->inceput;
        do
        {
            nr=nr+1;
            punct mij=mijloc_segment(m->p1, m->p2);
            X+=mij.x; Y+=mij.y; Z+=mij.z;
            m=m->urm;
        }while(m!=t->sfarsit->urm); (t->treiD).x=X/nr; (t->treiD).y=Y/nr; (t->treiD).z=Z/nr;

        if((t->treiD).x*eye.x+(t->treiD).y*eye.y+(t->treiD).z*eye.z<=D)
        {
            (t->doiD)=choose_pixel( intersectie_plan(t->treiD) );
            circle((t->doiD).x, (t->doiD).y, 3);
        }
        t=t->urm;
    }while(t!=NULL);
    setcolor(BLACK);
}

puncteCG* gasitCG(pereche s)
{
    puncteCG *t=primCG;
    if(primCG==NULL)return NULL;
    do
    {
        if(modul((t->doiD).x-s.x)<3 && modul((t->doiD).y-s.y)<3)return t;
        t=t->urm;
    }while(t!=NULL);return t;
}

void translare_corp(pereche p1, pereche p2, puncteCG* adresa)
{
    double difx=modul((double)(p1.x-p2.x));
    double dify=modul((double)(p1.y-p2.y)), L,K;

    punct p=(adresa->treiD);punct r;
    punct inter=intersectie_plan(p);

    double ei=dist_punct_punct(eye,inter), ip=dist_punct_punct(inter,p);

    if(p2.x>p1.x && p2.y<=p1.y)
    {
        punct p12=proiectie(inter,c1,c2);
        L=dify/(dist_punct_punct(inter,p12)-dify); K=(1.0)/(1+L);
        punct k={K*(inter.x+L*p12.x), K*(inter.y+L*p12.y), K*(inter.z+L*p12.z)};

        punct p23=proiectie(k,c2,c3);
        L=difx/(dist_punct_punct(k,p23)-difx); K=(1.0)/(1+L);
        r.x=K*(k.x+L*p23.x); r.y=K*(k.y+L*p23.y); r.z=K*(k.z+L*p23.z);
    }
    else if(p2.x<=p1.x && p2.y<p1.y)
    {
        punct p12=proiectie(inter,c1,c2);
        L=dify/(dist_punct_punct(inter,p12)-dify); K=(1.0)/(1+L);
        punct k={K*(inter.x+L*p12.x), K*(inter.y+L*p12.y), K*(inter.z+L*p12.z)};

        punct p14=proiectie(k,c1,c4);
        L=difx/(dist_punct_punct(k,p14)-difx); K=(1.0)/(1+L);
        r.x=K*(k.x+L*p14.x); r.y=K*(k.y+L*p14.y); r.z=K*(k.z+L*p14.z);
    }
    else if(p2.x<p1.x && p2.y>=p1.y)
    {
        punct p43=proiectie(inter,c4,c3);
        L=dify/(dist_punct_punct(inter,p43)-dify); K=(1.0)/(1+L);
        punct k={K*(inter.x+L*p43.x), K*(inter.y+L*p43.y), K*(inter.z+L*p43.z)};

        punct p14=proiectie(k,c1,c4);
        L=difx/(dist_punct_punct(k,p14)-difx); K=(1.0)/(1+L);
        r.x=K*(k.x+L*p14.x); r.y=K*(k.y+L*p14.y); r.z=K*(k.z+L*p14.z);
    }
    else
    {
        punct p43=proiectie(inter,c4,c3);
        L=dify/(dist_punct_punct(inter,p43)-dify); K=(1.0)/(1+L);
        punct k={K*(inter.x+L*p43.x), K*(inter.y+L*p43.y), K*(inter.z+L*p43.z)};

        punct p23=proiectie(k,c2,c3);
        L=difx/(dist_punct_punct(k,p23)-difx); K=(1.0)/(1+L);
        r.x=K*(k.x+L*p23.x); r.y=K*(k.y+L*p23.y); r.z=K*(k.z+L*p23.z);
    }
    double er=dist_punct_punct(eye,r), rpprim=er*ip/ei;
    L=er/rpprim; K=1+L;
    punct pprim={(r.x*K-eye.x)/L, (r.y*K-eye.y)/L, (r.z*K-eye.z)/L};

    double dif1=pprim.x-p.x, dif2=pprim.y-p.y, dif3=pprim.z-p.z;
    adresa->treiD=pprim;
    muchie*t=adresa->inceput;
    do
    {
        (t->p1).x+=dif1; (t->p2).x+=dif1; (t->p1).y+=dif2;
        (t->p2).y+=dif2; (t->p1).z+=dif3; (t->p2).z+=dif3;
        t=t->urm;
    }while(t!=adresa->sfarsit->urm);
}

void rotire_corp(pereche p1, pereche p2, puncteCG* &adr)
{
    muchie*t=adr->inceput;punct C=adr->treiD;
    do
    {
        (t->p1).x-=C.x; (t->p1).y-=C.y; (t->p1).z-=C.z;
        (t->p2).x-=C.x; (t->p2).y-=C.y; (t->p2).z-=C.z;
        t=t->urm;
    }while(t!=adr->sfarsit->urm);
    double L, K, distx, disty, U=0.000005; punct MIJ;

    if((p2.x>p1.x&&p1.y>=p2.y)||(p1.x>p2.x&&p2.y>=p1.y))///caz 1 si 4
    {
        if(p1.x>p2.x&&p2.y>=p1.y)U=-0.000005;

        distx=modul(p1.x-p2.x); disty=modul(p1.y-p2.y);
        L=disty/(Mpe2-disty);K=(1.0)/(1+L);
        punct k1={(centru_window.x+L*d2.x)*K, (centru_window.y+L*d2.y)*K, (centru_window.z+L*d2.z)*K};

        L=distx/(Npe2-distx);K=(1.0)/(1+L);
        punct k2={(centru_window.x+L*d3.x)*K, (centru_window.y+L*d3.y)*K, (centru_window.z+L*d3.z)*K};

        punct mij=mijloc_segment(k1,k2);
        plan P=ecuatie_plan(mij,centru_window,origine);

        ///refolosim mij -> devine normala pe plan

        mij.x=P.a; mij.y=P.b; mij.z=P.c; MIJ=mij;
        U*=sqrt((distx*distx+disty*disty)*(eye.x*eye.x+eye.y*eye.y+eye.z*eye.z));
    }
    else if((p1.x>=p2.x&&p1.y>p2.y)||(p2.x>=p1.x&&p2.y>p1.y))///caz 3 si 2
    {
        if(p2.x>=p1.x&&p2.y>p1.y)U=-0.000005;

        distx=modul(p1.x-p2.x); disty=modul(p1.y-p2.y);
        L=disty/(Mpe2-disty);K=(1.0)/(1+L);
        punct k1={(centru_window.x+L*d2.x)*K, (centru_window.y+L*d2.y)*K, (centru_window.z+L*d2.z)*K};

        L=distx/(Npe2-distx);K=(1.0)/(1+L);
        punct k2={(centru_window.x+L*d1.x)*K, (centru_window.y+L*d1.y)*K, (centru_window.z+L*d1.z)*K};

        punct mij=mijloc_segment(k1,k2);
        plan P=ecuatie_plan(mij,centru_window,origine);

        ///refolosim mij -> devine normala pe plan

        mij.x=P.a; mij.y=P.b; mij.z=P.c; MIJ=mij;
        U*=sqrt((distx*distx+disty*disty)*(eye.x*eye.x+eye.y*eye.y+eye.z*eye.z));
    }
    t=adr->inceput;
    do
    {
        t->p1=rotatie(t->p1,U,MIJ); (t->p1).x+=C.x; (t->p1).y+=C.y; (t->p1).z+=C.z;
        t->p2=rotatie(t->p2,U,MIJ); (t->p2).x+=C.x; (t->p2).y+=C.y; (t->p2).z+=C.z;
        t=t->urm;
    }while(t!=adr->sfarsit->urm);
}

void translare_punct(pereche p1, pereche p2, muchie* adr, punct &pun)
{
    double difx=modul((double)(p1.x-p2.x));
    double dify=modul((double)(p1.y-p2.y));double L,K;

    punct p=pun; punct r;
    punct inter=intersectie_plan(p);

    double ei=dist_punct_punct(eye,inter), ip=dist_punct_punct(inter,p);

    if(p2.x>p1.x && p2.y<=p1.y)
    {
        punct p12=proiectie(inter,c1,c2);
        L=dify/(dist_punct_punct(inter,p12)-dify); K=(1.0)/(1+L);
        punct k={K*(inter.x+L*p12.x), K*(inter.y+L*p12.y), K*(inter.z+L*p12.z)};

        punct p23=proiectie(k,c2,c3);
        L=difx/(dist_punct_punct(k,p23)-difx); K=(1.0)/(1+L);
        r.x=K*(k.x+L*p23.x); r.y=K*(k.y+L*p23.y); r.z=K*(k.z+L*p23.z);
    }
    else if(p2.x<=p1.x && p2.y<p1.y)
    {
        punct p12=proiectie(inter,c1,c2);
        L=dify/(dist_punct_punct(inter,p12)-dify); K=(1.0)/(1+L);
        punct k={K*(inter.x+L*p12.x), K*(inter.y+L*p12.y), K*(inter.z+L*p12.z)};

        punct p14=proiectie(k,c1,c4);
        L=difx/(dist_punct_punct(k,p14)-difx); K=(1.0)/(1+L);
        r.x=K*(k.x+L*p14.x); r.y=K*(k.y+L*p14.y); r.z=K*(k.z+L*p14.z);
    }
    else if(p2.x<p1.x && p2.y>=p1.y)
    {
        punct p43=proiectie(inter,c4,c3);
        L=dify/(dist_punct_punct(inter,p43)-dify); K=(1.0)/(1+L);
        punct k={K*(inter.x+L*p43.x), K*(inter.y+L*p43.y), K*(inter.z+L*p43.z)};

        punct p14=proiectie(k,c1,c4);
        L=difx/(dist_punct_punct(k,p14)-difx); K=(1.0)/(1+L);
        r.x=K*(k.x+L*p14.x); r.y=K*(k.y+L*p14.y); r.z=K*(k.z+L*p14.z);
    }
    else
    {
        punct p43=proiectie(inter,c4,c3);
        L=dify/(dist_punct_punct(inter,p43)-dify); K=(1.0)/(1+L);
        punct k={K*(inter.x+L*p43.x), K*(inter.y+L*p43.y), K*(inter.z+L*p43.z)};

        punct p23=proiectie(k,c2,c3);
        L=difx/(dist_punct_punct(k,p23)-difx); K=(1.0)/(1+L);
        r.x=K*(k.x+L*p23.x); r.y=K*(k.y+L*p23.y); r.z=K*(k.z+L*p23.z);
    }
    double er=dist_punct_punct(eye,r);
    double rpprim=er*ip/ei;

    L=er/rpprim; K=1+L;
    punct pprim={(r.x*K-eye.x)/L, (r.y*K-eye.y)/L, (r.z*K-eye.z)/L};

    muchie*t=adr;
    do
    {
        if((t->p1).x==pun.x && (t->p1).y==pun.y && (t->p1).z==pun.z)t->p1=pprim; else if((t->p2).x==pun.x && (t->p2).y==pun.y && (t->p2).z==pun.z)t->p2=pprim;
        t=t->urm;
    }while(t!=NULL);
    pun=pprim;
}

int char_la_int(char s[])
{
    int n=0, i=0, l=strlen(s);
    if(s[0] == '-')i++;
    for(int j=i;j<l;j++) n=n*10+s[j]-48;
    if(s[0]=='-')return -n; return n;
}

int numar_cifre(char s[])
{
    int nr_cif = 0;
    for(int i = 0; s[i]; i++)
        if(s[i] >= '0' && s[i] <= '9')nr_cif++;
    return nr_cif;
}

void citeste_sir(int x, int y, char mesaj[200], char S[200], bool isCifra)
{
    if(anulare) {return;}
    while (kbhit()) {
        getch();
    }
    setactivepage(1 - page);
    char Enter = 13, BackSpace = 8, Escape = 27, MultimeDeCaractereAcceptabile[200];
    if(isCifra) strcpy(MultimeDeCaractereAcceptabile, "-0123456789");
    else strcpy(MultimeDeCaractereAcceptabile, "-_0123456789abcdefghijklmnopqrstuvwxyz");
    char Sinitial[200], tasta;
    char Ss[200], mesajs[200];

    strcpy(Sinitial,S);
    settextjustify(LEFT_TEXT,TOP_TEXT);
    strcpy(mesajs,mesaj);
    strcat(mesajs,": ");
    if(light_selected)setbkcolor(fundal_intro);
    else setbkcolor(WHITE);
    outtextxy(x,y,mesajs);
    x=x+textwidth(mesajs);
    strcpy(S,"");
    do
    {
        tasta=getch();
        if (tasta==0)tasta=getch();
        else if (strchr(MultimeDeCaractereAcceptabile,tasta) && numar_cifre(Ss) <= 4)
        {
            strcpy(Ss,S);
            strcat(Ss,"_ ");
            setcolor(BLACK);
            outtextxy(x,y,Ss);
            char tt2[2];
            tt2[0]=tasta;
            tt2[1]=0;
            strcat(S,tt2);
            strcpy(Ss,S);
            strcat(Ss,"_ ");
            setcolor(BLACK);
            outtextxy(x,y,Ss);
        }
        else if (tasta == BackSpace)
        {
            if(strcmp(S,""))
            {
                strcpy(Ss,S);
                strcat(Ss,"_");
                setcolor(BLACK);
                outtextxy(x,y,Ss);
                setcolor(BLACK);
                S[strlen(S)-1]=0;
                strcpy(Ss,S);
                strcat(Ss,"_    ");
                outtextxy(x,y,Ss);
            }
        }
    }
    while (tasta!=Enter && tasta!=Escape);
    if (tasta == Escape) {anulare = true; tastare = false; return;}
    setcolor(BLACK);
    strcpy(Ss,S);
    strcat(Ss,"  ");
    outtextxy(x,y,Ss);
    setcolor(BLACK);
    outtextxy(x,y,S);
    setcolor(BLACK);
    tastare = false;
}

void back_to_hub()
{
    delay(100);
    alt_window_open2 = false;
    alt_window_open = false;
    closegraph(pagCrt);
    pagCrt = pagPrincipala;
    setcurrentwindow(pagPrincipala);
}

void open_fisier()
{
    char numeFis[256], numeFisAux[256], numeLabel[] = "Nume fisier";
    if(!romana_selected) strcpy(numeLabel, "File name");
    tastare = true;
    citeste_sir(60, 39, numeLabel, numeFisAux, 0);
    strcpy(numeFis, "saves/");
    strcat(numeFis, numeFisAux);
    strcat(numeFis, ".txt");
    if(anulare){anulare = false; back_to_hub(); return;}

    char semn[100];
    punct pct1; pct1.x = 0, pct1.y = 0, pct1.z = 0;
    ifstream fin(numeFis);
    fin.getline(semn, 100);
    if(!strcmp(semn, "eora5candfacasta")) adaugare_fisier_txt(pct1, 1, numeFis);
    else adaugare_fisier_obj(pct1, 100, numeFis);

    delay(100);
    fin.close();
    back_to_hub();
}

void save_fisier()
{
    char numeFis[256], numeFisAux[256], numeLabel[] = "Nume fisier";
    if(!romana_selected) strcpy(numeLabel, "File name");
    tastare = true;
    citeste_sir(60, 39, numeLabel, numeFisAux, 0);
    if(anulare){anulare = false; back_to_hub(); return;}

    strcpy(numeFis, "saves/");
    strcat(numeFis, numeFisAux);
    strcat(numeFis, ".txt");
    ofstream fout(numeFis);
    if(fout.is_open()) {
        fout<<"eora5candfacasta\n";
        puncteCG*g = primCG;
        while(g != NULL){
            muchie*t;
            t = new muchie; t = g->inceput;
            while(t != g->sfarsit->urm)
            {
                fout<<"m "<<t->p1.x<<" "<<t->p1.y<<" "<<t->p1.z<<" "<<t->p2.x<<" "<<t->p2.y<<" "<<t->p2.z<<"\n";
                t=t->urm;
            }
            fout<<"gc "<<g->treiD.x<<" "<<g->treiD.y<<" "<<g->treiD.z<<" "<<g->doiD.x<<" "<<g->doiD.y<<"\n";
            g = g -> urm;
        }
        fout.close();
    }
    else
        cout<<"eroare\n";

    back_to_hub();
}

void adaugare_fisier_txt(punct c, int L, char numeFisier[])
{
    char liniecrt[255];
    muchie*t, *I;
    puncteCG*g;
    ifstream fin(numeFisier);
    t = ultim;
    bool incep = true;
    while(fin.getline(liniecrt, sizeof(liniecrt)))
    {
        if(liniecrt[0] == 'm' && liniecrt[1] == ' ')
        {
            punct pct1, pct2;
            strcpy(liniecrt, liniecrt + 2);
            istringstream linie_citita(liniecrt);
            double numar;
            linie_citita>>numar; pct1.x = numar;
            linie_citita>>numar; pct1.y = numar;
            linie_citita>>numar; pct1.z = numar;
            linie_citita>>numar; pct2.x = numar;
            linie_citita>>numar; pct2.y = numar;
            linie_citita>>numar; pct2.z = numar;
            t=new muchie; ultim->urm=t; ultim=t; t->p1=pct1; t->p2=pct2; if(incep){I=t; incep=false;}//1
            ultim->urm=NULL;
        }
        else if(liniecrt[0] == 'g' && liniecrt[2] == ' ')
        {
            strcpy(liniecrt, liniecrt + 3);
            cout<<liniecrt<<"\n";
            istringstream linie_citita(liniecrt);
            double numar;
            if(primCG == NULL)
            {
                cout<<"gol\n";

                g=new puncteCG; primCG = g;
                ultimCG=primCG;
                g->inceput=I; g->sfarsit=ultim;
            }
            else
            {
                cout<<"nu e gol\n";
                g = new puncteCG; ultimCG->urm=g; ultimCG=g; g->inceput=I; g->sfarsit=ultim;
            }
            linie_citita>>numar; g->treiD.x = numar; linie_citita>>numar; g->treiD.y = numar;
            linie_citita>>numar; g->treiD.z = numar; linie_citita>>numar; g->doiD.x = numar;
            linie_citita>>numar; g->doiD.y = numar;

            incep = true;
            ultimCG->urm=NULL;
        }
    }
    fin.close();
}

void adaugare_fisier_obj(punct c, int L, char numeFisier[])
{
    char liniecrt[256];
    muchie*t, *I=ultim;
    punct pct1, pct2, pct3, pct4;int i = 1;
    ifstream fin(numeFisier);

    if(fin.getline(liniecrt, sizeof(liniecrt)))
        cout<<"merge\n";
    else
        cout<<"nu meerge!!!!!!!!!!!\n";
    while(fin.getline(liniecrt, sizeof(liniecrt)))
    {
        double numere[9999][4];
        if(liniecrt[0] == 'v' && liniecrt[1] == ' ')
        {
            int j = 1;
            char*ptr1 = liniecrt + 2;
            double nr = 0;
            cout<<liniecrt<<"\n";
            while(sscanf(ptr1, "%lf", &nr) == 1)
            {
                numere[i][j++] = nr;
                ptr1 = strchr(ptr1, ' ');
                if(!ptr1) break;
                ptr1++;
            }
            numere[i][1] *= L, numere[i][2] *= L, numere[i][3] *= L;
            numere[i][1] += c.x, numere[i][2] += c.y, numere[i][3] += c.z;
            i++;
        }
        else if(liniecrt[0] == 'f' && liniecrt[1] == ' ')
        {
            i = 0;
            char*ptr1 = liniecrt + 2;
            int nr = 0, bla[5], j = 1;
            while(sscanf(ptr1, "%d", &nr) == 1 && j <= 4)
            {
                cout<<nr<<" ";
                bla[j++] = nr;
                ptr1 = strchr(ptr1, ' ');
                if(!ptr1) break;
                ptr1++;
            }
            cout<<"\n";
            pct1.x = numere[bla[1]][1], pct1.y = numere[bla[1]][2], pct1.z = numere[bla[1]][3];
            pct2.x = numere[bla[2]][1], pct2.y = numere[bla[2]][2], pct2.z = numere[bla[2]][3];
            pct3.x = numere[bla[3]][1], pct3.y = numere[bla[3]][2], pct3.z = numere[bla[3]][3];
            if(j > 4){
                pct4.x = numere[bla[4]][1], pct4.y = numere[bla[4]][2], pct4.z = numere[bla[4]][3];
                t=new muchie; ultim->urm=t; ultim=t; t->p1=pct1; t->p2=pct2; //stage 1
                t=new muchie; ultim->urm=t; ultim=t; t->p1=pct2; t->p2=pct3; //stage 2
                t=new muchie; ultim->urm=t; ultim=t; t->p1=pct3; t->p2=pct4; // pentru 4
                t=new muchie; ultim->urm=t; ultim=t; t->p1=pct1; t->p2=pct4; // pentru 4
            }
            else{
                t=new muchie; ultim->urm=t; ultim=t; t->p1=pct1; t->p2=pct2; //stage 1
                t=new muchie; ultim->urm=t; ultim=t; t->p1=pct2; t->p2=pct3; //stage 2
            }
            ultim->urm=NULL; i++;
        }
    }

    double X=0,Y=0,Z=0,nrc=0;
    muchie*m=I->urm;
    do
    {
        punct mij=mijloc_segment(m->p1, m->p2);
        X+=mij.x; Y+=mij.y; Z+=mij.z;
        m=m->urm;
    }while(m!=NULL); punct cc={X/nrc, Y/nrc, Z/nrc};

    if(primCG==NULL){primCG=new puncteCG; primCG->treiD=c; primCG->doiD=choose_pixel(intersectie_plan(cc)); primCG->inceput=I->urm; primCG->sfarsit=ultim; ultimCG=primCG; ultimCG->urm=NULL;}
    else {ultimCG->urm=new puncteCG; ultimCG=ultimCG->urm; ultimCG->treiD=c; ultimCG->doiD=choose_pixel(intersectie_plan(cc)); ultimCG->inceput=I->urm; ultimCG->sfarsit=ultim; ultimCG->urm=NULL;}
}

void adaugare_linie(punct p1, punct p2)
{
    NRMUCHII++;
    punct centru_forma=mijloc_segment(p1,p2);
    muchie*I;muchie*t;
    t=new muchie; ultim->urm=t; ultim=t; t->p1=p1; t->p2=p2;I=t;
    ultim->urm=NULL;

    if(primCG==NULL){primCG=new puncteCG; primCG->treiD=centru_forma; primCG->doiD=choose_pixel(intersectie_plan(centru_forma)); primCG->inceput=I; primCG->sfarsit=ultim; ultimCG=primCG; ultimCG->urm=NULL;}
    else {ultimCG->urm=new puncteCG; ultimCG=ultimCG->urm; ultimCG->treiD=centru_forma; ultimCG->doiD=choose_pixel(intersectie_plan(centru_forma)); ultimCG->inceput=I; ultimCG->sfarsit=ultim; ultimCG->urm=NULL;}
}

void adaugare_paralelipiped(punct centru_forma, float lungime, float latime, float inaltime)
{
    muchie*t; punct k1=centru_forma, k2, k3, k4, k5, k6, k7, k8;
    muchie*I;
    NRMUCHII+=12;
    k1.y-=inaltime/2; k1.z+=latime/2; k1.x-=lungime/2;
    k2=k1; k2.z-=latime;
    k3=k2; k3.x+=lungime;
    k4=k3; k4.z+=latime;
    k5=k1; k5.y+=inaltime;
    k6=k5; k6.z-=latime;
    k7=k6; k7.x+=lungime;
    k8=k7; k8.z+=latime;

    t=new muchie; ultim->urm=t; ultim=t; t->p1=k1; t->p2=k4; I=t;//1
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k1; t->p2=k2;//2
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k4; t->p2=k3;//3
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k2; t->p2=k3;//4
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k1; t->p2=k5;//5
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k4; t->p2=k8;//6
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k5; t->p2=k8;//7
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k5; t->p2=k6;//8
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k2; t->p2=k6;//9
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k8; t->p2=k7;//10
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k3; t->p2=k7;//11
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k6; t->p2=k7;//12
    ultim->urm=NULL;

    if(primCG==NULL){primCG=new puncteCG; primCG->treiD=centru_forma; primCG->doiD=choose_pixel(intersectie_plan(centru_forma)); primCG->inceput=I; primCG->sfarsit=ultim; ultimCG=primCG; ultimCG->urm=NULL;}
    else {ultimCG->urm=new puncteCG; ultimCG=ultimCG->urm; ultimCG->treiD=centru_forma; ultimCG->doiD=choose_pixel(intersectie_plan(centru_forma)); ultimCG->inceput=I; ultimCG->sfarsit=ultim; ultimCG->urm=NULL;}
}

void adaugare_cilindru(punct c, float raza, int inaltime, int verticale)
{
    muchie*t;
    NRMUCHII+=(3*verticale);
    int i=0;
    muchie*I;
    for(i = 0; i < verticale; i++)
    {
        punct pct1, pct2, pct3;
        float tetha1 = (float)i * doipi / verticale;
        float tetha2 = (float)(i + 1) * doipi / verticale;

        pct1.x = raza * cos(tetha1);
        pct1.z = raza * sin(tetha1);
        pct1.y = 0;

        pct2.x = raza * cos(tetha2);
        pct2.z = raza * sin(tetha2);
        pct2.y = 0;

        pct3 = pct1; pct3.y = inaltime;

        pct1.x += c.x, pct1.y += c.y, pct1.z += c.z;
        pct2.x += c.x, pct2.y += c.y, pct2.z += c.z;
        pct3.x += c.x, pct3.y += c.y, pct3.z += c.z;

        punct pct1aux = pct1; pct1aux.y += inaltime;
        punct pct2aux = pct2; pct2aux.y += inaltime;

        t=new muchie; ultim->urm=t; ultim=t; t->p1=pct1; t->p2=pct2;if(i==0)I=t;
        t=new muchie; ultim->urm=t; ultim=t; t->p1=pct1; t->p2=pct3;
        t=new muchie; ultim->urm=t; ultim=t; t->p1=pct1aux; t->p2=pct2aux;
        ultim->urm=NULL;
    }

    double X=0,Y=0,Z=0,nr=0; muchie*m=I;
    do
    {
        nr=nr+1;
        punct mij=mijloc_segment(m->p1, m->p2);
        X+=mij.x; Y+=mij.y; Z+=mij.z;
        m=m->urm;
    }while(m!=NULL);

    punct rez={X/nr,Y/nr,Z/nr};
    double dif1=c.x-rez.x, dif2=c.y-rez.y, dif3=c.z-rez.z; m=I;
    do
    {
        (m->p1).x+=dif1; (m->p1).y+=dif2; (m->p1).z+=dif3;
        (m->p2).x+=dif1; (m->p2).y+=dif2; (m->p2).z+=dif3;
        m=m->urm;
    }while(m!=NULL);

    if(primCG==NULL){primCG=new puncteCG; primCG->treiD=c; primCG->doiD=choose_pixel(intersectie_plan(c)); primCG->inceput=I; primCG->sfarsit=ultim; ultimCG=primCG; ultimCG->urm=NULL;}
    else {ultimCG->urm=new puncteCG; ultimCG=ultimCG->urm; ultimCG->treiD=c; ultimCG->doiD=choose_pixel(intersectie_plan(c)); ultimCG->inceput=I; ultimCG->sfarsit=ultim; ultimCG->urm=NULL;}
}

void adaugare_con(punct c, float raza, int inaltime, int verticale)
{
    muchie*t;
    int i=0;
    muchie*I;
    NRMUCHII+=(2*verticale);
    for(i = 0; i < verticale; i++)
    {
        punct pct1, pct2, pct3;
        float tetha1 = (float)i * doipi / verticale;
        float tetha2 = (float)(i + 1) * doipi / verticale;

        pct1.x = raza * cos(tetha1);
        pct1.z = raza * sin(tetha1);
        pct1.y = 0;

        pct2.x = raza * cos(tetha2);
        pct2.z = raza * sin(tetha2);
        pct2.y = 0;

        pct3.x = 0, pct3.y = inaltime, pct3.z = 0;

        pct1.x += c.x, pct1.y += c.y, pct1.z +=  c.z;
        pct2.x += c.x, pct2.y += c.y, pct2.z +=  c.z;
        pct3.x += c.x, pct3.y += c.y, pct3.z +=  c.z;

        t=new muchie; ultim->urm=t; ultim=t; t->p1=pct1; t->p2=pct2;if(i==0)I=t;
        t=new muchie; ultim->urm=t; ultim=t; t->p1=pct1; t->p2=pct3;
        ultim->urm=NULL;
    }

    double X=0,Y=0,Z=0,nr=0; muchie*m=I;
    do
    {
        nr=nr+1;
        punct mij=mijloc_segment(m->p1, m->p2);
        X+=mij.x; Y+=mij.y; Z+=mij.z;
        m=m->urm;
    }while(m!=NULL);

    punct rez={X/nr,Y/nr,Z/nr};
    double dif1=c.x-rez.x, dif2=c.y-rez.y, dif3=c.z-rez.z; m=I;
    do
    {
        (m->p1).x+=dif1; (m->p1).y+=dif2; (m->p1).z+=dif3;
        (m->p2).x+=dif1; (m->p2).y+=dif2; (m->p2).z+=dif3;
        m=m->urm;
    }while(m!=NULL);

    if(primCG==NULL){primCG=new puncteCG; primCG->treiD=c; primCG->doiD=choose_pixel(intersectie_plan(c)); primCG->inceput=I; primCG->sfarsit=ultim; primCG->urm=NULL; ultimCG=primCG;}
    else {ultimCG->urm=new puncteCG; ultimCG=ultimCG->urm; ultimCG->treiD=c; ultimCG->doiD=choose_pixel(intersectie_plan(c)); ultimCG->inceput=I; ultimCG->sfarsit=ultim; ultimCG->urm=NULL;}
}

void adaugare_sfera(punct c, float raza, int verticale, int orizontale)
{
    muchie*t,*I;
    NRMUCHII+=(2*verticale*orizontale);
    for(int i = 0; i < verticale; i++)
    {
        float tetha1 = (float)i * doipi / verticale;
        float tetha2 = (float)(i + 1) * doipi / verticale;

        float cost1=cos(tetha1);
        float cost2=cos(tetha2);
        float sint1=sin(tetha1);
        float sint2=sin(tetha2);

        float unusupraorizontale=(1.0)/orizontale;

        for(int j = 0; j < orizontale; j++)
        {
            punct pct1, pct2, pct3;
            float phi1 = (float)j * pi * unusupraorizontale;
            float phi2 = (float)(j + 1) * pi *unusupraorizontale;

            float cosf1=cos(phi1);
            float cosf2=cos(phi2);
            float sinf1=sin(phi1);
            float sinf2=sin(phi2);

            pct1.x = raza * cost2 * sinf1;
            pct1.y = raza * sint2 * sinf1;
            pct1.z = raza * cosf1;

            pct2.x = raza * cost2 * sinf2;
            pct2.y = raza * sint2 * sinf2;
            pct2.z = raza * cos(phi2);

            pct3.x = raza * cost1 * sinf2;
            pct3.y = raza * sint1 * sinf2;
            pct3.z = raza * cosf2;

            pct1.x += c.x, pct1.y += c.y, pct1.z +=  c.z;
            pct2.x += c.x, pct2.y += c.y, pct2.z +=  c.z;
            pct3.x += c.x, pct3.y += c.y, pct3.z +=  c.z;

            t=new muchie; ultim->urm=t; ultim=t; t->p1=pct1; t->p2=pct2;if(i==0&&j==0)I=t;
            t=new muchie; ultim->urm=t; ultim=t; t->p1=pct2; t->p2=pct3;

            ultim->urm=NULL;
        }
    }

    double X=0,Y=0,Z=0,nr=0; muchie*m=I;
    do
    {
        nr=nr+1;
        punct mij=mijloc_segment(m->p1, m->p2);
        X+=mij.x; Y+=mij.y; Z+=mij.z;
        m=m->urm;
    }while(m!=NULL);

    punct rez={X/nr,Y/nr,Z/nr};
    double dif1=c.x-rez.x, dif2=c.y-rez.y, dif3=c.z-rez.z; m=I;
    do
    {
        (m->p1).x+=dif1; (m->p1).y+=dif2; (m->p1).z+=dif3;
        (m->p2).x+=dif1; (m->p2).y+=dif2; (m->p2).z+=dif3;
        m=m->urm;
    }while(m!=NULL);

    if(primCG==NULL){primCG=new puncteCG; primCG->treiD=c; primCG->doiD=choose_pixel(intersectie_plan(c)); primCG->inceput=I; primCG->sfarsit=ultim; ultimCG=primCG; ultimCG->urm=NULL;}
    else {ultimCG->urm=new puncteCG; ultimCG=ultimCG->urm; ultimCG->treiD=c; ultimCG->doiD=choose_pixel(intersectie_plan(c)); ultimCG->inceput=I; ultimCG->sfarsit=ultim; ultimCG->urm=NULL;}
}

void adaugare_cub(punct centru_forma, float size_cub)
{
    NRMUCHII+=12;
    float latura=size_cub;
    size_cub /= 2;
    muchie*I;
    muchie*t; punct k1=centru_forma, k2, k3, k4, k5, k6, k7, k8;

    k1.y-=size_cub; k1.z+=size_cub; k1.x-=size_cub;
    k2=k1; k2.z-=latura;
    k3=k2; k3.x+=latura;
    k4=k3; k4.z+=latura;
    k5=k1; k5.y+=latura;
    k6=k5; k6.z-=latura;
    k7=k6; k7.x+=latura;
    k8=k7; k8.z+=latura;

    t=new muchie; ultim->urm=t; ultim=t; t->p1=k1; t->p2=k4; I=t;//1
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k1; t->p2=k2;     //2
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k4; t->p2=k3;     //3
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k2; t->p2=k3;     //4
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k1; t->p2=k5;     //5
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k4; t->p2=k8;     //6
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k5; t->p2=k8;     //7
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k5; t->p2=k6;     //8
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k2; t->p2=k6;     //9
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k8; t->p2=k7;     //10
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k3; t->p2=k7;     //11
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k6; t->p2=k7;     //12
    ultim->urm=NULL;

    if(primCG==NULL){primCG=new puncteCG; primCG->treiD=centru_forma; primCG->doiD=choose_pixel(intersectie_plan(centru_forma)); primCG->inceput=I; primCG->sfarsit=ultim; ultimCG=primCG; ultimCG->urm=NULL;}
    else {ultimCG->urm=new puncteCG; ultimCG=ultimCG->urm; ultimCG->treiD=centru_forma; ultimCG->doiD=choose_pixel(intersectie_plan(centru_forma)); ultimCG->inceput=I; ultimCG->sfarsit=ultim; ultimCG->urm=NULL;}
}

void adaugare_tetraedru(punct p, float latura)
{
    NRMUCHII+=6;
    muchie*I;
    punct k1, k2, k3, k4, baza;
    double hpe2=latura/sqrt(6);
    k4=p; k4.y+=hpe2;
    k2=p; k2.y-=hpe2; baza=k2; k1=k2; double b=latura/sqrt(3); k2.z+=b;
    k1.z-=(b*0.5); k1.x-=(latura*0.5); k3=k1; k3.x+=latura;

    muchie*t;
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k1; t->p2=k2; I=t;
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k2; t->p2=k3;
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k3; t->p2=k1;
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k4; t->p2=k1;
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k4; t->p2=k2;
    t=new muchie; ultim->urm=t; ultim=t; t->p1=k4; t->p2=k3;
    ///t=new muchie; ultim->urm=t; ultim=t; t->p1=k4; t->p2=baza;
    ultim->urm=NULL;

    double X=0,Y=0,Z=0,nr=0; muchie*m=I;
    do
    {
        nr=nr+1;
        punct mij=mijloc_segment(m->p1, m->p2);
        X+=mij.x; Y+=mij.y; Z+=mij.z;
        m=m->urm;
    }while(m!=NULL);

    punct rez={X/nr,Y/nr,Z/nr};
    double dif1=p.x-rez.x, dif2=p.y-rez.y, dif3=p.z-rez.z; m=I;
    do
    {
        (m->p1).x+=dif1; (m->p1).y+=dif2; (m->p1).z+=dif3;
        (m->p2).x+=dif1; (m->p2).y+=dif2; (m->p2).z+=dif3;
        m=m->urm;
    }while(m!=NULL);

    if(primCG==NULL){primCG=new puncteCG; primCG->treiD=p; primCG->doiD=choose_pixel(intersectie_plan(p)); primCG->inceput=I; primCG->sfarsit=ultim; ultimCG=primCG; ultimCG->urm=NULL;}
    else {ultimCG->urm=new puncteCG; ultimCG=ultimCG->urm; ultimCG->treiD=p; ultimCG->doiD=choose_pixel(intersectie_plan(p)); ultimCG->inceput=I; ultimCG->sfarsit=ultim; ultimCG->urm=NULL;}
}

void coordonate_forme_pres(int shape)
{
    char sX[10],sY[10],sZ[10],sL[10],x[]="x",y[]="y",z[]="z",l[15], latim[15];
    char sX1[10],sY1[10],sZ1[10],x1[]="x1",y1[]="y1",z1[]="z1";
    char sX2[10],sY2[10],sZ2[10],x2[]="x2",y2[]="y2",z2[]="z2";
    char sR[10],sVert[10],sOriz[10],r[15],vert[15],oriz[15],inalt[15];
    char sInalt[10], sLatime[10];
    tastare = true;
    if(!romana_selected) {strcpy(l,"Length");strcpy(latim,"Width");strcpy(r,"Radius");strcpy(vert,"Strips");strcpy(oriz, "Layers");strcpy(inalt, "Height");}
    else {strcpy(l,"Lungime");strcpy(latim,"Latime");strcpy(r,"Raza");strcpy(vert,"Fasii");strcpy(oriz, "Straturi");strcpy(inalt, "Inaltime");}

    switch(shape)
    {
    case 1:
        citeste_sir(160, 13, x, sX, 1); citeste_sir(160, 42, y, sY, 1); citeste_sir(160, 72, z, sZ, 1);
        citeste_sir(160, 104, l, sL, 1);
        break;
    case 2:
        citeste_sir(160, 13, x, sX, 1); citeste_sir(160, 42, y, sY, 1); citeste_sir(160, 72, z, sZ, 1);
        citeste_sir(160, 104, l, sL, 1);
        break;
    case 3:
        citeste_sir(160, 13, x, sX, 1); citeste_sir(160, 42, y, sY, 1); citeste_sir(160, 72, z, sZ, 1);
        citeste_sir(160, 104, r, sR, 1); citeste_sir(160, 133, vert, sVert, 1); citeste_sir(160, 163, oriz, sOriz, 1);
        break;
    case 4:
    case 5:
        citeste_sir(160, 13, x, sX, 1); citeste_sir(160, 42, y, sY, 1); citeste_sir(160, 72, z, sZ, 1);
        citeste_sir(160, 104, r, sR, 1); citeste_sir(160, 133, inalt, sInalt, 1); citeste_sir(160, 163, vert, sVert, 1);
        break;
    case 6:
        citeste_sir(160, 13, x, sX, 1); citeste_sir(160, 42, y, sY, 1); citeste_sir(160, 72, z, sZ, 1);
        citeste_sir(160, 104, l, sL, 1); citeste_sir(160, 133, latim, sLatime, 1); citeste_sir(160, 163, inalt, sInalt, 1);
        break;
    case 7:
        citeste_sir(160, 13, x1, sX1, 1); citeste_sir(160, 42, y1, sY1, 1); citeste_sir(160, 72, z1, sZ1, 1);
        citeste_sir(160, 104, x2, sX2, 1); citeste_sir(160, 133, y2, sY2, 1); citeste_sir(160, 163, z2, sZ2, 1);
    default:
        break;
    }

    tastare = false;
    alt_window_open = false;
    closegraph(pagForme);
    pagCrt = pagPrincipala;
    setcurrentwindow(pagPrincipala);

    if(anulare){anulare = false; return;}

    punct centru; centru.x = 0; centru.y = 0; centru.z = 0;
    punct pct1, pct2;
    pct1.x = char_la_int(sX1), pct1.y = char_la_int(sY1), pct1.z = char_la_int(sZ1);
    pct2.x = char_la_int(sX2), pct2.y = char_la_int(sY2), pct2.z = char_la_int(sZ2);
    centru.x = char_la_int(sX);
    centru.y = char_la_int(sY);
    centru.z = char_la_int(sZ);
    float L = char_la_int(sL);
    float R = char_la_int(sR);
    float H = char_la_int(sInalt);
    int fasii =  char_la_int(sVert);
    int straturi = char_la_int(sOriz);
    int latime = char_la_int(sLatime);

    switch(shape)
    {
    case 1:
        adaugare_cub(centru, L);
        break;
    case 2:
        adaugare_tetraedru(centru, L);
        break;
    case 3:
        adaugare_sfera(centru, R, fasii, straturi);
        break;
    case 4:
        adaugare_con(centru, R, H, fasii);
        break;
    case 5:
        adaugare_cilindru(centru, R, H, fasii);
        break;
    case 6:
        adaugare_paralelipiped(centru, L, latime, H);
        break;
    case 7:
        adaugare_linie(pct1, pct2);
    default:
        break;
    }
}

void check_input_overlay()
{
    int x, y;
    if(kbhit())
    {
        int tasta = getch();
        if(tasta == 27)
        {
            cout<<"se inchide "<<pagCrt<<"\n";
            alt_window_open = false;
            if(pagCrt == pagPrincipala) exit(0);
            back_to_hub();
        }
    }
    if(ismouseclick(WM_LBUTTONDOWN))
    {
        getmouseclick(WM_LBUTTONDOWN, x, y);
        clearmouseclick(WM_LBUTTONDOWN);
        if(intro)
        {
            if(x >= 200 && x <= 300 && y >= 300 && y <= 350){
                intro = false;
            }
            else if(x >= 50 && x <= 125 && y >= 60 && y <= 110){
                romana_selected = true;
            }
            else if(x >= 135 && x <= 210 && y >= 60 && y <= 110){
                romana_selected= false;
            }
            else if(x >= 290 && x <= 365 && y >= 60 && y <= 110){
                light_selected = true;
            }
            else if(x >= 375 && x <= 450 && y >= 60 && y <= 110){
                light_selected = false;
            }
            else if(x >= 165 && x <= 215 && y >= 200 && y <= 250){
                color = 0;
            }
            else if(x >= 225 && x <= 275 && y >= 200 && y <= 250){
                color = 1;
            }
            else if(x >= 285 && x <= 335 && y >= 200 && y <= 250){
                color = 2;
            }

        }
        else {
            if(alt_window_open && !alt_window_open2) {
                if(x >= 10 && x <= 150 && y >= 10 && y <= 40){
                    coordonate_forme_pres(1); delay(100);
                }
                else if(x >= 10 && x <= 150 && y >=40 && y <= 70){
                    coordonate_forme_pres(2); delay(100);
                }
                else if(x >= 10 && x <= 150 && y >= 70 && y <= 100) {
                    coordonate_forme_pres(3); delay(100);
                }
                else if(x >= 10 && x <= 150 && y >= 100 && y <= 130) {
                    coordonate_forme_pres(4); delay(100);
                }
                else if(x >= 10 && x <= 150 && y >= 130 && y <= 160) {
                    coordonate_forme_pres(5); delay(100);
                }
                else if(x >= 10 && x <= 150 && y >= 160 && y <= 190) {
                    coordonate_forme_pres(6); delay(100);
                }
                else if(x >= 10 && x <= 150 && y >= 190 && y <= 220) {
                    coordonate_forme_pres(7); delay(100);
                }
            }
            else if(!alt_window_open && !alt_window_open2)
            {
                if(x >= 0 && x <= 100 && y >= 0 && y <= 30){
                    char save[10];
                    if(romana_selected) strcpy(save, "Salvare");
                    else strcpy(save, "Save");
                    pagSave = initwindow(600, 100, save, 400, 300), setcurrentwindow(pagSave), pagCrt = pagSave;
                    settextstyle(3, 0, 2); alt_window_open = true, alt_window_open2 = true;
                    for(int i = 1; i <= 3; i++) {rutina_test(1), page = 1 - page;} rutina_test(1);
                    delay(100), save_fisier();
                    }
                else if(x >= 100 && x <= 200 && y >= 0 && y <= 30){
                    char open[10];
                    if(romana_selected) strcpy(open, "Deschide");
                    else strcpy(open, "Open");
                    pagOpen = initwindow(600, 100, open, 400, 300), setcurrentwindow(pagOpen), pagCrt = pagOpen;
                    settextstyle(3, 0, 2), alt_window_open = true, alt_window_open2 = true;
                    for(int i = 1; i <= 3; i++) {rutina_test(1), page = 1 - page;} rutina_test(1);
                    delay(50), open_fisier();
                }
                else if(x >= 200 && x <= 300 && y >= 0 && y <= 30){
                    stergere();
                }
                else if(x >= 300 && x <= 400 && y >= 0 && y <= 30){
                    char forme[10];
                    if(romana_selected) strcpy(forme, "Forme");
                    else strcpy(forme, "Shapes");
                    pagForme = initwindow(320, 230, forme, 400, 200); setcurrentwindow(pagForme), pagCrt = pagForme;
                    settextstyle(3, 0, 2); alt_window_open = true;
                }
            }
        }
    }
}

void gui_save_open()
{
    if(light_selected) {
        culoare_butoane = fundal_intro; setfillstyle(SOLID_FILL, fundal_intro);}
    else{
        culoare_butoane = WHITE; setfillstyle(SOLID_FILL, WHITE);}
    rectangle(50, 30, 550, 70);
    floodfill(51, 31, BLACK);
    setbkcolor(fundal);
}

void gui_forme_menu()
{
    char cub[10]="Cub", tetraedru[15]="Tetraedru", sfera[10]="Sfera", con[10]="Con";
    char cilindru[10]="Cilindru",  paralalelipiped[25]="Paralalelipiped", linie[10]="Linie";

    if(light_selected) {
        culoare_butoane = fundal_intro; setfillstyle(SOLID_FILL, fundal_intro);}
    else{
        culoare_butoane = WHITE; setfillstyle(SOLID_FILL, WHITE);}

    rectangle(10, 10, 310, 40);
    floodfill(11, 11, BLACK);

    rectangle(10, 40, 310, 70);
    floodfill(11, 41, BLACK);

    rectangle(10, 70, 310, 100);
    floodfill(11, 71, BLACK);

    rectangle(10, 100, 310, 130);
    floodfill(11, 101, BLACK);

    rectangle(10, 130, 310, 160);
    floodfill(11, 131, BLACK);

    rectangle(10, 160, 310, 190);
    floodfill(11, 161, BLACK);

    rectangle(10, 190, 310, 220);
    floodfill(11, 191, BLACK);

    if(!romana_selected) {
        strcpy(cub,"Cube");strcpy(tetraedru,"Tetrahedron");strcpy(sfera,"Sphere");
        strcpy(con,"Cone");strcpy(cilindru,"Cylinder");strcpy(paralalelipiped, "Parallelepiped");
        strcpy(linie,"Line");
        setbkcolor(culoare_butoane);
        outtextxy(20, 14, cub); outtextxy(20, 43, tetraedru);outtextxy(20, 73, sfera);
        outtextxy(20, 103, con); outtextxy(20, 133, cilindru);outtextxy(20, 162, paralalelipiped);
        outtextxy(20, 194, linie);
    }
    else {
        setbkcolor(culoare_butoane);
        outtextxy(20, 14, cub); outtextxy(20, 43, tetraedru);outtextxy(20, 73, sfera);
        outtextxy(20, 103, con); outtextxy(20, 133, cilindru);outtextxy(20, 163 , paralalelipiped);
        outtextxy(20, 194, linie);
    }
    line(150, 10, 150, 220);
    setbkcolor(fundal);
    if(!tastare) check_input_overlay();
}

void gui_menus()
{
    char save[10]="Salvare", open[10]="Deschide", reset[10]="Resetare", forme[10]="Forme";

    if(light_selected) {
        culoare_butoane = fundal_intro; setfillstyle(SOLID_FILL, fundal_intro);}
    else{
        culoare_butoane = WHITE; setfillstyle(SOLID_FILL, WHITE);}

    rectangle(0, 0, 100, 30);
    floodfill(1, 1, BLACK);
    rectangle(100, 0, 200, 30);
    floodfill(101, 1, BLACK);
    rectangle(200, 0, 300, 30);
    floodfill(201, 1, BLACK);
    rectangle(300, 0, 400, 30);
    floodfill(301, 1, BLACK);

    if(!romana_selected) {
        strcpy(save,"Save"); strcpy(open,"Open");strcpy(reset,"Reset"); strcpy(forme,"Shapes ");
        setbkcolor(culoare_butoane);
        outtextxy(30, 4, save); outtextxy(126, 4, open); outtextxy(226, 4, reset); outtextxy(320, 4, forme);
    }
    else {
        setbkcolor(culoare_butoane);
        outtextxy(20, 4, save); outtextxy(112, 4, open); outtextxy(215, 4, reset); outtextxy(322, 4, forme);
    }
    setbkcolor(fundal);
    if(!tastare) check_input_overlay();
}

void gui_intro()
{
    char start[20]="START", limba[20]="Limba:", tema[20]="Tema:", culoare[30]="Culoare linii:";
    if(!romana_selected) {
        strcpy(limba,"Language:"); strcpy(tema,"Theme:");strcpy(culoare,"Line color:");}
    setfillstyle(SOLID_FILL, WHITE);
    rectangle(200, 300, 300, 350);
    floodfill(201, 301, BLACK);
    setbkcolor(WHITE);
    outtextxy(220, 314, start);
    setbkcolor(fundal_intro);
    setcolor(BLACK);

    outtextxy(105, 30, limba);
    setcolor(WHITE);
    rectangle(50, 60, 125, 110);
    readimagefile("romana.jpg",51,61,124,109);
    rectangle(135, 60, 210, 110);
    readimagefile("engleza.jpg",136,61,209,109);
    setcolor(WHITE);
    if(romana_selected){for(int i = 0; i < 3; i++) rectangle(47-i, 57-i, 128+i, 113+i);}
    else {for(int i = 0; i < 3; i++) rectangle(132-i, 57-i, 213+i, 113+i);}

    setcolor(BLACK);
    outtextxy(345, 30, tema);
    setfillstyle(SOLID_FILL, COLOR(204, 230, 252));
    rectangle(290, 60, 365, 110);
    floodfill(291, 61, BLACK);
    setfillstyle(SOLID_FILL, COLOR(47,49,54));
    rectangle(375, 60, 450, 110);
    floodfill(376, 61, BLACK);
    setcolor(WHITE);
    if(light_selected){for(int i = 0; i < 3; i++) rectangle(287-i, 57-i, 368+i, 113+i); fundal = COLOR(204, 230, 252);}
    else {for(int i = 0; i < 3; i++) rectangle(372-i, 57-i, 453+i, 113+i); fundal = COLOR(47,49,54);}

    setcolor(BLACK);
    outtextxy(200, 170, culoare);
    setfillstyle(SOLID_FILL, BLACK);
    rectangle(165, 200, 215, 250);
    floodfill(166, 201, BLACK);
    setfillstyle(SOLID_FILL, COLOR(255,0,0));
    rectangle(225, 200, 275, 250);
    floodfill(226, 201, BLACK);
    setfillstyle(SOLID_FILL, COLOR(255,192,203));
    rectangle(285, 200, 335, 250);
    floodfill(286, 201, BLACK);
    setcolor(WHITE);
    if(color == 0){for(int i = 0; i < 3; i++) rectangle(162-i, 197-i, 218+i, 253+i); culoare_linie = BLACK;}
    else if(color == 1){for(int i = 0; i < 3; i++) rectangle(222-i, 197-i, 278+i, 253+i); culoare_linie = COLOR(255,0,0);}
    else if(color == 2){for(int i = 0; i < 3; i++) rectangle(282-i, 197-i, 338+i, 253+i); culoare_linie = COLOR(255,192,203);}

    if(!tastare)check_input_overlay();
}

void rutina_test(bool start_to_project_from_beggining)
{
    ///0- incep de la prim->urm->urm->urm
    ///1-incep normal, cu axele alea rosu, verge si albastru
    setactivepage(page);
    setvisualpage(1-page);
    cleardevice();
    if(!alt_window_open){if(!intro) gui_menus();else gui_intro();}
    else {if(!alt_window_open2) gui_forme_menu();else gui_save_open();}
    delay(1);
    proiectare(start_to_project_from_beggining);
    proiectareCG();
}

int main()
{
    initializare();
    pagIntro = initwindow(500, 400, "Intro");
    settextstyle(3, 0, 2);
    setbkcolor(fundal);
    gui_intro();
    while(intro)
    {rutina_test(1);page = 1-page;}

    initializare();
    pagPrincipala = initwindow(N, M, "PROIECT 3D");
    pagCrt = pagPrincipala;
    settextstyle(3, 0, 2);
    proiectare(1);
    setbkcolor(fundal);
    gui_menus();

    POINT cursorpos;
    puncteCG*A;
    while(true)
    {
        GetCursorPos(&cursorpos);
        rutina_test(1);
        page = 1 - page;
        pereche sageata={mousex(), mousey()};
        if(ismouseclick(WM_LBUTTONDOWN) && !tastare)
        {
            sageata.x=mousex(); sageata.y=mousey();
            A=gasitCG(sageata);
            if(A!=NULL)
            {
                GetCursorPos(&cursorpos);
                while(GetAsyncKeyState(VK_LBUTTON))
                {
                    rutina_test(1);
                    POINT P; GetCursorPos(&P);
                    if(P.x==cursorpos.x && P.y==cursorpos.y && GetAsyncKeyState(VK_MBUTTON)){ }
                    pereche p1={cursorpos.x,cursorpos.y} ,p2={P.x,P.y};
                    translare_corp(p1,p2,A);
                    cursorpos=P;
                    page=1-page;
                }
            }
            else
            {
                pereche R1, R2;
                punct PUN={0,0,0};
                muchie*t=prim->urm->urm->urm;
                if(t!=NULL)
                do
                {
                    R1=choose_pixel(intersectie_plan(t->p1));
                    R2=choose_pixel(intersectie_plan(t->p2));
                    if((modul(R1.x-sageata.x)<4 && modul(R1.y-sageata.y<4)) || (modul(R2.x-sageata.x)<4 && modul(R2.y-sageata.y<4)))break;
                    t=t->urm;

                }while(t!=NULL);

                if(t!=NULL)
                {
                    GetCursorPos(&cursorpos);
                    if(modul(R1.x-sageata.x)<2 && modul(R1.y==sageata.y)<2) PUN=t->p1;else PUN=t->p2;
                    if(modul(choose_pixel(intersectie_plan(PUN)).x-sageata.x)<4 && modul(choose_pixel(intersectie_plan(PUN)).y-sageata.y)<4)
                    while(GetAsyncKeyState(VK_LBUTTON))
                    {
                        rutina_test(1);
                        POINT P; GetCursorPos(&P);
                        pereche p1={cursorpos.x,cursorpos.y} ,p2={P.x,P.y};
                        translare_punct(p1,p2,t,PUN);
                        cursorpos=P;
                        page=1-page;
                    }
                }
            }
        }
        else if(ismouseclick(WM_MBUTTONDOWN) && !tastare)
        {
            sageata.x=mousex(); sageata.y=mousey();
            A=gasitCG(sageata);
            if(A!=NULL)
            {
                GetCursorPos(&cursorpos);
                while(GetAsyncKeyState(VK_MBUTTON))
                {
                    rutina_test(1);
                    POINT P; GetCursorPos(&P);
                    pereche p1={cursorpos.x,cursorpos.y} ,p2={P.x,P.y};
                    rotire_corp(p1,p2,A);
                    cursorpos=P;
                    page=1-page;
                }
            }
        }
        else if(GetAsyncKeyState(VK_BACK) && !tastare)
        {
            sageata.x=mousex(); sageata.y=mousey();
            A=gasitCG(sageata);
            if(A!=NULL)
            {
                muchie*t=prim;
                do
                {
                    t=t->urm;
                }while(t->urm!=A->inceput);t=A->sfarsit;

                if(A==primCG&&primCG->urm==NULL){primCG=ultimCG=NULL;}
                //else if(A==primCG&&primCG!=ultimCG){primCG=primCG->urm;}
                /*else
                {
                    puncteCG*u=primCG;
                    do
                    {
                        u=u->urm;
                    }while(u->urm!=A);u->urm=A->urm;
                }*/
                rutina_test(1);
                page=1-page;
            }
        }

        while(GetAsyncKeyState(VK_RBUTTON) && !tastare)
        {
            rutina_test(1);
            POINT P; GetCursorPos(&P);
            pereche p1={cursorpos.x,cursorpos.y} ,p2={P.x,P.y};
            ajustare_puncte_principale_la_rotatie(p1,p2);
            cursorpos=P;
            page = 1 - page;
        }

        while(GetAsyncKeyState(VK_UP) && dist_punct_punct(eye,origine)>22.2 && !tastare)
        {
            rutina_test(1);
            ajustare_puncte_principale_la_zoom(22.2,1);
            page = 1 - page;
        }

        while(GetAsyncKeyState(VK_DOWN) && !tastare)
        {
            rutina_test(1);
            ajustare_puncte_principale_la_zoom(22.2,0);
            page = 1-page;
        }

        while(GetAsyncKeyState(VK_RIGHT) && !tastare)
        {
            rutina_test(1);

            eye=rotatie_y(eye);
            centru_window=rotatie_y(centru_window);

            c1=rotatie_y(c1);
            c2=rotatie_y(c2);
            c3=rotatie_y(c3);
            c4=rotatie_y(c4);

            d1=mijloc_segment(c1,c4);
            d2=mijloc_segment(c1,c2);
            d3=mijloc_segment(c2,c3);
            d4=mijloc_segment(c3,c4);

            //window=ecuatie_plan(c1,c2,c3);
            page = 1 - page;
        }

        while(GetAsyncKeyState(VK_CONTROL) && !tastare)
        {

        rutina_test(0);
        page=1-page;

        }
    }
    getch();closegraph();
    return 0;
}
