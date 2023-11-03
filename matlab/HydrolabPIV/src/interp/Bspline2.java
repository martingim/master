import java.lang.Thread;
import java.io.Serializable;

public class Bspline2 implements Serializable {
    private double[] tx;
    private double[] ty;
    private double[][] c;
    private int I;
    private int J;
    private int k;

    private static int maxthreads = Integer.MAX_VALUE;

    transient private double[][] xi;
    transient private double[][] yi;
    transient private double[][] f;
    transient private double[][] idx;

    public Bspline2(int n) {
	I = 4;
	J = 4;
	c = new double[I][J];
	tx = new double[I+4];
	ty = new double[I+4];
	for(int i=I; i < I+4; i++) {
	    tx[i]=1.0;
	    ty[i]=1.0;
	}
	k = 4;
	maxthreads = n;
    }
    public Bspline2(double[][] coeff, double[] tx,double[] ty, int order) {
	this.c = coeff;
	this.tx = tx;
	this.ty = ty;
	this.I = coeff.length;
	this.J = coeff[0].length;
	this.k = order;
    }
    

    public static void setMaxThreads(int n) {
	maxthreads = n;
    }

    public static double[][] basis(int i, int k, double[][] x, double[] t) {
	
	int I = t.length-k-1;//+1;
	int M = x.length;
	int N = x[0].length;
	double[][] rs = new double[M][N];

	for(int m=0; m<M; m++) {
	    for(int n=0; n<N; n++) {
		rs[m][n] = basis(i-1,k,x[m][n],t,I);
	    }
	}
	return rs;
    }

    private static double basis(int i, int k, double x, double[] t,int I) {

	// Outside interval
	if (i<0) return 0.0;
	if (i>I) return 0.0;

	if ((i==I) && essentiallyEqual(x,t[I+k-1])) return 1.0;       

	if (k==1) {
	    if ((t[i] <= x) && (x < t[i+1]))
		return 1.0;
	    else
		return 0.0;
	} else {
	    double N = 0.0;
	    if (!essentiallyEqual(t[i+k-1],t[i])) {
		N = (x-t[i])/(t[i+k-1]-t[i])*basis(i,k-1,x,t,I);      
	    }    
	    if (!essentiallyEqual(t[i+k],t[i+1])) {       
		N = N + (t[i+k]-x)/(t[i+k] - t[i+1])*basis(i+1,k-1,x,t,I);
	    }    
	    return N;
	}
	
    }
   
    private static boolean essentiallyEqual(double a, double b) {       
	double epsilon = 1e-10;
	if (Math.abs(a) > Math.abs(b))
	    return Math.abs(a - b) <= Math.abs(b) * epsilon;
	else	    
	    return Math.abs(a - b) <= Math.abs(a) * epsilon;
    }
 
    private class Worker implements Runnable {
	private int s;
	private int S;
			
	public Worker(int s, int S) {
	    this.s = s;
	    this.S = S;
	}

	public void run() {
	    int M = xi.length;
	    int N = xi[0].length;

	    int m,n,i,j,i0,j0,m0,n0,np,mp;;
	    double x,y;

	    if(M<N) {
		m0 = 0;
		mp = 1;
		n0 = s;
		np = S;
	    } else {
		m0 = s;
		mp = S;
		n0 = 0;
		np = 1;		
	    }

	    for(m = m0; m < M; m += mp) {    
		for(n = n0; n < N; n += np) {
		    x = xi[m][n];
		    y = yi[m][n];
		    
		    if ((tx[0] <= x) && (x <= tx[J+k-1]) && (ty[0] <= y) && (y <= ty[I+k-1])) {
			idx[m][n] = 1.0;

			// There is probably a faster way to do this
			i0 = k-1; 
			for(i = k; i < I; i++)
			    if (y >= ty[i]) i0++;
			j0 = k-1; 
			for(j = k; j < J; j++)
			    if (x >= tx[j]) j0++;		
			
			for(i = -k+1; i < 1; i++) {
			    for(j = -k+1; j < 1; j++) {	    			   
				f[m][n] += c[i0+i][j0+j]*basis(i0+i,k,yi[m][n],ty,I-1)*basis(j0+j,k,xi[m][n],tx,J-1);
			    }	  	  
			}      
		    }
		    //f[m][n] = s; debug
		}
	    }
	    
	}
    }
    
    public double[][] mask() {
	return idx;
    }

    public double[][] evaluate(double[][] xi, double[][] yi) {
	//int I = t.length-k-1;//+1;
	int M = xi.length;
	int N = xi[0].length;
	f = new double[M][N];
	idx = new double[M][N];
	this.xi = xi;
	this.yi = yi;
	int S = Runtime.getRuntime().availableProcessors();
	S = (S<=maxthreads) ? S : maxthreads;

	if(S==1) {
	    int m,n,i,j,i0,j0;
	    double x,y;

	    for(m = 0; m < M; m++) {    
		for(n = 0; n < N; n++) {
		    x = xi[m][n];
		    y = yi[m][n];
		    
		    if ((tx[0] <= x) && (x <= tx[J+k-1]) && (ty[0] <= y) && (y <= ty[I+k-1])) {
			idx[m][n] = 1.0;

			// There is probably a faster way to do this
			i0 = k-1; 
			for(i = k; i < I; i++)
			    if (y >= ty[i]) i0++;
			j0 = k-1; 
			for(j = k; j < J; j++)
			    if (x >= tx[j]) j0++;		
			
			for(i = -k+1; i < 1; i++) {
			    for(j = -k+1; j < 1; j++) {	    			   
				f[m][n] += c[i0+i][j0+j]*basis(i0+i,k,yi[m][n],ty,I-1)*basis(j0+j,k,xi[m][n],tx,J-1);
			    }	  	  
			}      
		    }
		    //f[m][n] = s; debug
		}
	    }
	} else {

	    Thread[] thread = new Thread[S];
	    for(int s=0; s<S; s++) {
		thread[s] = new Thread(new Worker(s,S));
		thread[s].start();
	    }
	    for(int s=0; s<S; s++) {	
		try {  
		    thread[s].join();
		} catch (java.lang.InterruptedException e) {
		}
	    }
	}
	return f;
    }
}
