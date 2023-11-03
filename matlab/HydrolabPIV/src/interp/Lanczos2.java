import java.lang.Thread;

public class Lanczos2 {
    private double[][] c;
    private int I;
    private int J;
    transient private double[][] x;
    transient private double[][] y;
    transient private double[][] f;
    transient private double[][] mask;
    transient private int a;

    private static int maxthreads = Integer.MAX_VALUE;

    public static void setMaxThreads(int n) {
	maxthreads = n;
    }

    public Lanczos2(int n) {
	I=1;
	J=1;
	c = new double[I][J];
	maxthreads = n;	
    }

    public Lanczos2(double[][] coeff) {
	c = coeff;
	I = coeff.length;
	J = coeff[0].length;
    }
    
    private double lfun(double x,double a) {	
	if(x==0) {
	    return 1.0;
	} else if(Math.abs(x)>a) {
	    return 0.0;
	} else {
	    return a*Math.sin(Math.PI*x)*Math.sin(Math.PI*x/a)/(Math.PI*Math.PI*x*x);
	}	
    }

    public double[][] evaluate(double[][] x, double[][] y, int a) {
	
	int M = x.length;
	int N = x[0].length;

	f = new double[M][N];
	mask = new double[M][N];
	
	int S = Runtime.getRuntime().availableProcessors();
	S = (S<=maxthreads) ? S : maxthreads;

	if(S==1) {
	    int m,n,i,j;
	    double xi,yi;       

	    for(m = 0; m < M; m++) {
		for(n = 0; n < N; n++) {
		    if(1 > y[m][n] || y[m][n] > I || 1 > x[m][n] || x[m][n] > J) {
			f[m][n] = 0.0;
			mask[m][n] = 0.0;
		    } else {		    
			xi = Math.floor(x[m][n]) - 1;
			yi = Math.floor(y[m][n]) - 1;
			for(i=(int)yi-a+1; i < (int)yi+a+1; i++) {
			    for(j=(int)xi-a+1; j < (int)xi+a+1; j++) {
				if (i < 0) continue;
				if (j < 0) continue;
				if (i >= I) continue;
				if (j >= J) continue;
				f[m][n] += c[i][j]*lfun(y[m][n]-i-1,(double) a)*lfun(x[m][n]-j-1,(double) a);
			    }   
			}
			mask[m][n] = 1.0;
		    }
		}
	    }
	} else {
	    this.x = x;
	    this.y = y;
	    this.a = a;

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
    
    

	
    public double[][] mask() {
	return this.mask;
    }
    
    private class Worker implements Runnable {
	private int s;
	private int S;
			
	public Worker(int s, int S) {
	    this.s = s;
	    this.S = S;
	}

	public void run() {
	    int M = x.length;
	    int N = x[0].length;

	    int m,n,i,j,m0,n0,np,mp;
	    double xi,yi;       
	    
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
	    	    if(1 > y[m][n] || y[m][n] > I || 1 > x[m][n] || x[m][n] > J) {
			f[m][n] = 0.0;
			mask[m][n] = 0.0;
		    } else {		    
			xi = Math.floor(x[m][n]) - 1;
			yi = Math.floor(y[m][n]) - 1;
			for(i=(int)yi-a+1; i < (int)yi+a+1; i++) {
			    for(j=(int)xi-a+1; j < (int)xi+a+1; j++) {
				if (i < 0) continue;
				if (j < 0) continue;
				if (i >= I) continue;
				if (j >= J) continue;
				f[m][n] += c[i][j]*lfun(y[m][n]-i-1,(double) a)*lfun(x[m][n]-j-1,(double) a);
			    }   
			}
			mask[m][n] = 1.0;
		    }
		}
	    }
	}
    }
    
}
