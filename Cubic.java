import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;

import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;


public class Cubic {

    public static void progonka(double[] A, double[] B, double[] C, double[] F, int N, double[] c, double[] f, double[] h, double epsilon, int iteration){
        double[]alpha = new double[N+1];
        double[]betta = new double[N+1];

        alpha[2]= B[1]/C[1];//14110.;//B[1]/C[1];//0.;//-C[1]/F[1];//0.;//B[0]/C[0];
        betta[2]= F[1]/C[1];//25.;//F[1]/C[1];//0.;//h/F[1];//0.;//f[0]/C[0];
        for(int i=2;i<N;i++)
        {
            alpha[i+1]=B[i]/(C[i]-alpha[i]*A[i]);
            betta[i+1]=(A[i]*betta[i]+F[i])/(C[i]-alpha[i]*A[i]);
        }
        alpha[1]=(C[1]*alpha[2]-B[1])/(alpha[2]*A[1]);
        betta[1]=(betta[2]*(C[1]-alpha[1]*A[1]) - F[1])/A[1];

//        alpha[1] = B[1]/C[1];
//        betta[1] = F[1]/C[1];
//        for(int i=2;i<N;i++)
//        {
//            alpha[i]=B[i]/(C[i]+alpha[i-1]*A[i]);
//            betta[i]=(F[i]-A[i]*betta[i-1])/(C[i]+alpha[i-1]*A[i]);
//        }
//        A[N] = h[N];
//        C[N] = (-2.)*h[N];
//        F[N] = (f[N]-f[N-1])/h[N];
//        alpha[N] = 0.;
//        betta[N] = (F[N]-A[N]*betta[N-1])/(C[N]+alpha[N-1]*A[N]);



//        System.out.println(alpha[N]);
//        System.out.println(betta[N]);

//        System.out.println([N]);

//        0 вариант:
//        c[0]=2.;//-Math.PI*Math.PI;
//        c[N]=2.;//Math.PI*Math.PI;
//        1 вариант:

//        c[0] = 0.;
//        c[N] = 30.;
//        c[0]=0.; //вторая производная sin(Pi*x/2) в 0
//        c[N]=(-Math.PI*Math.PI)/4.;//вторая производная sin(Pi*x/2) в 1
//        c[N]=-Math.PI*Math.PI*Math.sqrt(2.)/32.; //вторая производная cos(Pi*x/4) в 0
//        c[0]=(-Math.PI*Math.PI)/16.;//вторая производная cos(Pi*x/4) в 1
//        c[0]=-Math.PI*Math.PI/4.+1./(epsilon*epsilon); //вторая производная cos(Pi*x/4)+epsilon в 0
//        c[N]=Math.exp(-1./epsilon)/(epsilon*epsilon);//вторая производная cos(Pi*x/4)+epsilon в 1
//        c[N]=-Math.cos(1.5); //вторая производная cos(Pi*x/4) в 0
//        c[0]=-Math.cos(0.5);//вторая производная cos(Pi*x/4) в 1
        c[0]=0.; //вторая производная sinx в 0
        c[N]=-Math.sin(1.);//вторая производная sinx в 1

//        c[0] = 0.;
//        c[N] = 2.*Math.tan(1.)/(Math.cos(1.)*Math.cos(1.));

//        c[0] = 0.;
//        c[N] = 5.;
//          System.out.println("Вторая производная Math.sin(Math.PI*x/2.) + Math.exp(-x/epsilon)");
//          c[0]=1./(epsilon*epsilon);//вторая производная Math.sin(Math.PI*x/2.) + Math.exp(-x/epsilon) в 0
//          c[N]=(-Math.PI*Math.PI)/4. + Math.exp(-1./epsilon)/(epsilon*epsilon);//вторая производная Math.sin(Math.PI*x/2.) + Math.exp(-x/epsilon) в 1

//        2 вариант:
//        double proizv1 = Math.PI/2.;//производная sin(Pi*x/2) в 0
//        double proizv2 = 0.;//производная sin(Pi*x/2) в 1
//        double proizv1 = 0.;//производная cos(Pi*x/2) в 0
//        double proizv2 = -Math.PI/2.;//производная cos(Pi*x/2) в 1
//        double proizv1 = 1.;//производная sin(x) в 0
//        double proizv2 = Math.cos(1.);//производная sin(x) в 1
//        double proizv1 = 0.;//производная cos(x) в 0
//        double proizv2 = -Math.sin(1.);//производная cos(x) в 1
//        double proizv1 = 0.;//производная sin(x)*sin(x) в 0
//        double proizv2 = 2.*Math.sin(1.)*Math.cos(1.);//производная sin(x)*sin(x) в 1


//        System.out.println("Первая производная Math.sin(Math.PI*x/2.) + Math.exp(-x/epsilon)");
//        double proizv1 = Math.PI/2.-1./epsilon; //производная Math.sin(Math.PI*x/2.) + Math.exp(-x/epsilon) в 0
//        double proizv2 = - (Math.exp(-1./epsilon)/epsilon); //производная Math.sin(Math.PI*x/2.) + Math.exp(-x/epsilon) в 1

//          double proizv1 = -1./epsilon; //производная Math.cos(Math.PI*x/2.) + Math.exp(-x/epsilon) в 0
//          double proizv2 = -Math.PI/2.- (Math.exp(-1.0/epsilon)/epsilon); //производная Math.cos(Math.PI*x/2.) + Math.exp(-x/epsilon) в 1

//        double proizv1 = -1./epsilon;
//        double proizv2 = -Math.PI/2. - Math.exp(-1./epsilon)/epsilon;
//        c[0] = 6.*alpha[1]/(h[1]*(2.*alpha[1]+1.)) * ((h[1]*betta[1])/(6.*alpha[1]) + (f[1]-f[0])/h[1] - proizv1);
//        c[N] = 6.*((f[N-1]-f[N])/h[N] - h[N]*betta[N]/6. - proizv2)/(h[N]*(2.+alpha[N]));

//        System.out.println(c[N]);

//        iteration = 0;
        for (int i=N-1;i>0;i--)
        {
            c[i]=alpha[i+1]*c[i+1] + betta[i+1];
//            System.out.println(c[i]);
            iteration++;
        }

//        System.out.println("iteration = "+iteration);
    }


    public static  double[] ravnomSetka(int N){
        double[] h = new double[N+1];
        for (int i=1;i<N+1;i++){
            h[i] = 1./N;
        }
        return h;
    }

    //    public static void setkaShishkina(double epsilon, int N, double[] uzelx/*, double[] uzely*/, double sigma1, /*double sigma2,*/ double[] h/*, double[] hk*/){
    public static double[] setkaShishkina(double epsilon, int N){
//        double ShishkinMesh(double epsilon, int N, double *uzelx, double *uzely, double sigma1, double sigma2, double *h, double *hk){
        double sigma1;
        double[] h = new double[N+1];
        sigma1 = Math.min(0.5, 4.*epsilon*Math.log(N*1.0)/1.);
//            sigma2 = Math.min(0.5, 2.*epsilon*Math.log(N*1.0)/2.);

        for(int i=1;i<=N;i++){
            if ((i<=N/2)&&(i>=1)) {
                h[i] = 2.*sigma1/N;
//                    hk[i] = 2.*sigma2/N;
            }
            if ((i<=N)&&(i>N/2)) {
                h[i] = 2.*(1.-sigma1)/N;
//                    hk[i] = 2.*(1.-sigma2)/N;
            }
//                cout<<"h["<<i<<"] = "<<h[i]<<endl;
//                cout<<"hk["<<i<<"] = "<<hk[i]<<endl;
        }

        double uzelx[] = new double[N+1];
        uzelx[0]=0.;
//            uzely[0]=0.;
        for (int i=1;i<N+1;i++)
        {
            uzelx[i]=uzelx[i-1]+h[i];
//                uzely[i]=uzely[i-1]+hk[i];

        }

//            for(int i=1;i<=N;i++)
//            {
//                cout<<"uzelx"<<i<<" = "<<uzelx[i]<<endl;
//                cout<<"uzely"<<i<<" = "<<uzely[i]<<endl;
//            }
//            return 0;
        return h;

    }

    public static double[] setkaBakhvalova(double epsilon, int N){
        double uzelx[] = new double[N+1];
        double sigma1;
        double[] h = new double[N+1];
        sigma1 = Math.min(0.5,(-4.)*epsilon*Math.log(epsilon)/1.);
//        sigma2 = min(0.5,(-2.)*epsilon*log(epsilon)/2.);
        if (epsilon>Math.exp(-1)) {
            sigma1 = 0.5;
//            sigma2 = 0.5;
        }


        if (sigma1==0.5){
//            cout<<"NULL FIVE!"<<endl;
            uzelx[0]=0.;
            for(int i=1;i<=N;i++){
                h[i]=1./N;
                uzelx[i]=uzelx[i-1]+h[i];
            }
        }
//        if (sigma2==0.5){
//            uzely[0]=0.;
//            for(int i=1;i<=N;i++){
//                hk[i]=1./N;
//                uzely[i]=uzely[i-1]+hk[i];
//            }
//        }

        if (sigma1<0.5)
        {
            uzelx[0] = 0.;
            for (int i=1;i<=N/2;i++)
            {
                uzelx[i] = (-2.)*epsilon*Math.log(1.-2.*(1.-epsilon)*i/N);
                //cout<<uzelx[i]<<endl;
            }
            for (int i=N/2;i<=N;i++)
            {
                uzelx[i] = sigma1 + (2.*i/N - 1.)*(1.-sigma1);
            }
            for(int i=1;i<=N;i++)
            {
                h[i]=uzelx[i]-uzelx[i-1];
            }
        }

//        if (sigma2<0.5)
//        {
//            uzely[0] = 0.;
//            for (int i=1;i<=N/2;i++)
//            {
//                uzely[i] = (-2.)*epsilon*log(1-2*(1-epsilon)*i/N)/2.;
//                //cout<<uzely[i]<<endl;
//            }
//            for (int i=N/2.;i<=N;i++)
//            {
//                //cout<<"i = "<<i*2./N<<endl;
//                uzely[i] = sigma2 + (2.*i/N - 1.)*(1.-sigma2);
//                cout<<"uzely = "<<uzely[i]<<endl;
//            }
//            for(int i=1;i<=N;i++)
//            {
//                hk[i]=uzely[i]-uzely[i-1];
//                //cout<<"hy = "<<hk[i]<<endl;
//            }
//        }
        return h;
    }

    public static void findPoints(int N, double[] h, double[] uzel){
        uzel[0]= 0.;
        for (int i=1;i<N+1;i++)
        {
            uzel[i] = uzel[i-1]+h[i];
//            System.out.println(uzel[i]);
        }
    }


    public static void findFunction(int N, double[] f, double[] uzel, Function<Double, Double> function){
        for (int i=0;i<N+1;i++)
        {
            f[i] = function.apply(uzel[i]);
//            System.out.println(i+" = "+f[i]);
        }
    }

    public static void findCoeff(double[] h, double[] f, double[] A, double[] B, double[] C, double[] F, int N){
        for (int i = 1; i < N; i++) {
            A[i] = h[i];
            C[i] = (-2.)*(h[i] + h[i+1]);
            B[i] = h[i+1];
            F[i] = (-6.) * ( (f[i + 1] - f[i])/h[i+1] - (f[i] - f[i - 1])/h[i] );
//            F[i] = (-3.) * ( (f[i + 1] - f[i])/h[i+1] - (f[i] - f[i - 1])/h[i] );
//            F[i] = -6. * (f[i + 1] - 2. * f[i] + f[i - 1]) / (h);
//            F[i] = 3. * (f[i + 1] - 2. * f[i] + f[i - 1]) / (h);
//            F[i] = f[i];
        }
//        F[0] = -6.*(f[1] - 2*f[0]) / h;
//        B[0] = h;
//        C[0] = (-4.)*h;
//        A[0] = h;

    }








    public static void findS(int L, int N, double[] uzel, double[] x, double[] f, double[] S, double[] h, double[] hh, double[] c){
//        double hh=1./L;
//        double[] hh = ravnomSetka(L);
        x[0]=0.;

        for(int k=1;k<L+1;k++){
            x[k]=x[k-1]+hh[k];
        }

        for(int i=0;i<N+1;i++){
            for (int k=0;k<L;k++){
                // if ( ((x[k]>=uzel[i])&&(x[k]<uzel[i+1]))||(x[k]==1.) )
                if  ((x[k]>=uzel[i])&&(x[k]<uzel[i+1]))
                {

                    //  cout<<"x["<<k<<"] ="<<x[k]<<" popal v "<<"[ uzel["<<i<<"]; uzel["<<i+1<<"] ) = "<<uzel[i]<<";"<<uzel[i+1]<<endl;
                    //cout<<x[k]<<" popal v "<<"[" <<uzel[i]<<";"<< uzel[i+1]<< ")"<<endl;
                    //if (k==L) {cout<<"!!!!"<<endl;}
//                    S[k]=  f[i]+ ( (f[i]-f[i-1])/h + h*(2*c[i]+c[i-1])/3. )*(x[k]-uzel[i]) + c[i] * (x[k]-uzel[i])*(x[k]-uzel[i]) + ( (c[i]-c[i-1])/(3.*h) )*(x[k]-uzel[i])*(x[k]-uzel[i])*(x[k]-uzel[i]);
//                    S[k]=  f[i+1]+ ( (f[i+1]-f[i])/h + h*(2*c[i+1]+c[i])/3 )*(x[k]-uzel[i]) + c[i+1] * (x[k]-uzel[i])*(x[k]-uzel[i]) + ( (c[i+1]-c[i])/(3*h) )*(x[k]-uzel[i])*(x[k]-uzel[i])*(x[k]-uzel[i]);

//                    S[k] = (uzel[i+1] - x[k])*(uzel[i+1] - x[k])*(uzel[i+1] - x[k])*c[i]/(6*h) + (x[k] - uzel[i])*(x[k] - uzel[i])*(x[k] - uzel[i])*c[i+1]/(6*h) + (f[i+1]/h - c[i+1]*h/6.)*(x[k] - uzel[i]) + (f[i]/h - c[i]*h/6.)*(uzel[i+1] - x[k]);

                    S[k] = (uzel[i+1] - x[k])*(uzel[i+1] - x[k])*(uzel[i+1] - x[k])*c[i]/(6.*h[i+1]) + (x[k] - uzel[i])*(x[k] - uzel[i])*(x[k] - uzel[i])*c[i+1]/(6.*h[i+1]) + (f[i+1]/h[i+1] - c[i+1]*h[i+1]/6.)*(x[k] - uzel[i]) + (f[i]/h[i+1] - c[i]*h[i+1]/6.)*(uzel[i+1] - x[k]);
//
//                    System.out.println("S["+k+"] = "+S[k]);
                }
            }
        }

        S[L] = (uzel[N] - x[L])*(uzel[N] - x[L])*(uzel[N] - x[L])*c[N]/(6.*h[N]) + (x[L] - uzel[N-1])*(x[L] - uzel[N-1])*(x[L] - uzel[N-1])*c[N]/(6.*h[N]) + (f[N]/h[N] - c[N]*h[N]/6.)*(x[L] - uzel[N-1]) + (f[N-1]/h[N] - c[N-1]*h[N]/6.)*(uzel[N] - x[L]);

//        SplineInterpolator interpolator = new SplineInterpolator();
//        PolynomialSplineFunction spline = interpolator.interpolate(uzel, f);
    }

    public static void pogreshnost(double[] res, int L, double[] x, double[] S, int N, Function<Double, Double> function, List<Double> norm, double[] uzel, double[] f){
        AkimaSplineInterpolator interpolator = new AkimaSplineInterpolator();
        PolynomialSplineFunction spline = interpolator.interpolate(uzel, f);
        double norma = 0.;
        for(int i=0;i<L;i++)
//        for (int i=L/N;i<=L-L/N;i++)
//            for(int i=N;i<=L-N;i++)
        {
//            res[i] = Math.abs(function.apply(x[i]) - S[i]);
            res[i] = Math.abs(function.apply(x[i]) - spline.value(x[i]));
//            System.out.println("res["+i+"] = "+res[i]);
        }

        //норма погрешности:
        for (int i=0;i<L+1;i++)
//        for(int i=N;i<=L-L/N;i++)
//        for (int i=L/N;i<=L-L/N;i++)
//        for(int i=N;i<=L-N;i++)
        {
            if(res[i]>norma) norma=res[i];
//            System.out.println(i+" = "+norma);
        }
        norm.add(norma);
//        iteration++;
//        System.out.println("norma "+norma);
//        cout<<"norma = "<<norma<<endl;// вывод нормы погрешности
        //cout<<"epsilon = "<<epsilon<<endl<<" N = "<<N<<endl;
        //cout<<"iteration: "<<iteration<<endl;
//        cout<<"N = "<<N<<endl;

//        System.out.println("N = "+N);
    }

//    public static void memory(){
//
//    }


    public static double[] cubicalSpline(int N, Function<Double, Double> function, double epsilon, List<Double> norm){
        final int L = (5*N);
        double[] x = new double[L+1];

        double[]A = new double[N+1];
        double[]B = new double[N+1];
        double[]C = new double[N+1];
        double[]F = new double[N+1];
        double[]f = new double[N+1];
        double[]uzel = new double[N+1];
        double[]c = new double[N+1];

        double[]res = new double[L+1];
        double[]S = new double[L+1];
        double[] hRavnom = ravnomSetka(N);
        double[] hRavnomGust = ravnomSetka(L);
//        double h = 1./N;
//        double epsilon = 10.e-2;
        double[] hSetkaShishkina = setkaShishkina(epsilon, N);
        double[] hSetkaShishkinaGust = setkaShishkina(epsilon, L);
        double[] hSetkaBakhvalova = setkaBakhvalova(epsilon, N);
        double[] hSetkaBakhvalovaGust = setkaBakhvalova(epsilon, L);

        int iteration = 0;
//        double[]h = hSetkaBakhvalova;
//        double[]hGust = hSetkaBakhvalovaGust;
        double[]h = hSetkaShishkina;
        double[]hGust = hSetkaShishkinaGust;
//        double[]h = hRavnom;
//        double[]hGust = hRavnomGust;

        findPoints(N, h, uzel);
        findFunction(N, f, uzel, function);
        findCoeff(h, f, A, B, C, F, N);
        progonka(A, B, C, F, N, c, f, h, epsilon, iteration);
        findS(L, N, uzel, x, f, S, h, hGust, c);
        pogreshnost(res, L, x, S, N, function, norm, uzel, f);
//        System.out.println("iteration = "+iteration);
        return S;
    }



    public static void bicubicalSpline(int N, BiFunction<Double, Double, Double> biFunction){
        double[] Y = new double[N+1];
        double[]uzel = new double[N+1];
        double[]h = new double[N+1];
        double epsilon = 1.e-2;
        h = ravnomSetka(N);
        findPoints(N,h, uzel);
        List<Double> norm = new ArrayList<>();
        for(int i=0;i<N+1;i++){
            double finalI = Double.valueOf(i);
            Function<Double, Double> function = x -> biFunction.apply(x, finalI);
            cubicalSpline(N, function, epsilon, norm);
        }
        System.out.println("Максимальная норма погрешности " + Collections.max(norm));
        System.out.println("Минимальная норма погрешности " + Collections.min(norm));
    }

    public static void main(String[] args){
//        System.out.println("Enter the number of nodes:");
//        Scanner s = new Scanner(System.in);
//        int N = s.nextInt();
//        List<double> norm = new ArrayList<>();
//        bicubicalSpline(100, (x,y)->Math.sin(Math.PI*x/2.) + Math.exp(y));
//        cubicalSpline(N, x->Math.sin(Math.PI*x/2.), norm);
//        System.out.println("norma = "+norm);
        Map<Integer, List> map = new HashMap<>();
        double epsilon = 1.e-5;
        System.out.println("Epsilon = "+epsilon);
        double a = 0.;
        for (int i=16; i<=4096; i=i*2){
            List<Double> norm = new ArrayList<>();
//            cubicalSpline(i, x->x*x*x*x*x, epsilon, norm);
//            cubicalSpline(i, x->Math.cos(x+0.5), epsilon, norm);
//            cubicalSpline(i, x->Math.cos(Math.PI*x/4.), epsilon, norm);
//            cubicalSpline(i, x->Math.sin(Math.PI*x/2.) + Math.exp(-x/epsilon), epsilon, norm);
//            cubicalSpline(i, x->Math.sin(x)*Math.sin(x), epsilon, norm);
//            cubicalSpline(i, x->Math.sin(x), epsilon, norm);
            cubicalSpline(i, x->Math.cos(Math.PI*x/2.) + Math.exp(-x/epsilon), epsilon, norm);
//            cubicalSpline(i, x->Math.sin(Math.PI*x/2.) + Math.exp(-x/epsilon), epsilon, norm);
            double four = a/norm.get(0);
            a = norm.get(0);
//            System.out.println("||"+i+"|| = "+ norm);
            System.out.println("Порядок точности log2 (||"+i/2.+"||/||"+i+"||) = "+Math.log10(four)/Math.log10(2.));
            System.out.println("_______________________________________________________________");
            System.out.println("||"+i+"|| = "+ norm);

//            map.put(i, norm);
        }
//        Collections.sort(map.keySet());
//        System.out.println(map);
    }
}

