//import java.lang.Thread;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

public class DiffMeasure 
{
    private double[][] f1;
    private double[][] m1;
    private double[][] f2;
    private double[][] m2;
    private int[] idh;
    private int[] idw;
    private double alpha;
    private double[][] cc;
    private double[][] mm;
    private boolean normalized;
    
    private static ExecutorService executor;
    private static int nThreads;
    

    public DiffMeasure() 
    {
	if(executor==null) {
	    nThreads =  Runtime.getRuntime().availableProcessors();
	    executor = Executors.newFixedThreadPool(nThreads);
	}

    }
    
    public DiffMeasure(int n) 
    {        	
	nThreads = n;
	executor = Executors.newFixedThreadPool(n);
    }

    public double[][] maskedncc(double[][] f1, double[][] m1, double[][] f2, double[][] m2, int[] idh, int[] idw) 
    {
        return maskednmd(f1,m1,f2,m2,idh,idw,Double.NaN);
    }

    public double[][] maskedncc(double[][] f1, double[][] m1, double[][] f2, double[][] m2, int[] idh, int[] idw, double alpha) 
    {        
        return maskednmd(f1,m1,f2,m2,idh,idw,Double.NaN);
    }
    
    public double[][] maskedcc(double[][] f1, double[][] m1, double[][] f2, double[][] m2, int[] idh, int[] idw) 
    {
        return maskedmd(f1,m1,f2,m2,idh,idw,Double.NaN);
    }

    public double[][] maskedcc(double[][] f1, double[][] m1, double[][] f2, double[][] m2, int[] idh, int[] idw, double alpha) 
    {        
        return maskedmd(f1,m1,f2,m2,idh,idw,Double.NaN);
    }
    
    public double[][] maskedmd(double[][] f1, double[][] m1, double[][] f2, double[][] m2, int[] idh, int[] idw, double alpha) 
    {
        int ma = f1.length;
        int na = f1[0].length;
        int mb = f1.length;
        int nb = f1[0].length;        

        this.cc = new double[ma+mb-1][na+nb-1];
        this.mm = new double[ma+mb-1][na+nb-1];
	
	this.f1 = f1;
	this.m1 = m1;
	this.f2 = f2;
	this.m2 = m2;
	this.idh = idh;
	this.idw = idw;
	this.alpha = alpha;       

	ExecutorCompletionService<Integer> ecs = new ExecutorCompletionService<Integer>(executor);

	List<Future<Integer>> list = new ArrayList<Future<Integer>>(); 
 
	for(int s=0; s<nThreads; s++) {	    
	    Callable<Integer> callable = new Worker1(s,nThreads);
	    Future<Integer> future  = ecs.submit(callable);	   	  	   
	    list.add(future);
	}
	for(Future<Integer> future : list){	
	    try {
		Integer result = future.get();			    	    		
	    } catch (InterruptedException | ExecutionException e) {
	    }
	}	
        
        return cc;     
    }
    
    public double[][] maskednmd(double[][] f1, double[][] m1, double[][] f2, double[][] m2, int[] idh, int[] idw, double alpha) 
    {	
        int ma = f1.length;
        int na = f1[0].length;
        int mb = f1.length;
        int nb = f1[0].length;
        
        this.cc = new double[ma+mb-1][na+nb-1];
        this.mm = new double[ma+mb-1][na+nb-1];
		
    	
	this.f1 = f1;
	this.m1 = m1;
	this.f2 = f2;
	this.m2 = m2;
	this.idh = idh;
	this.idw = idw;
	this.alpha = alpha;
    
	ExecutorCompletionService<Integer> ecs = new ExecutorCompletionService<Integer>(executor);

	List<Future<Integer>> list = new ArrayList<Future<Integer>>(); 
 
	for(int s=0; s<nThreads; s++) {	    
	    Callable<Integer> callable = new Worker2(s,nThreads);
	    Future<Integer> future  = ecs.submit(callable);	   	  	   
	    list.add(future);
	}
	for(Future<Integer> future : list){	
	    try {
		Integer result = future.get();			    	    		
	    } catch (InterruptedException | ExecutionException e) {
	    }
	}
    
        return cc;
    }
       
    private class Worker1 implements Callable<Integer>
    {
        private int s;
        private int S;
	
        public Worker1(int s, int S) 
        { 
    	    this.s = s;
            this.S = S;
    	}

        public Integer call() throws Exception
        {
            int ma = f1.length;
            int na = f1[0].length;
            int mb = f1.length;
            int nb = f1[0].length;

    	    int I = idh.length;
            int J = idw.length;
            int i0,ip,j0,jp;

            if(I<J) {
                i0 = 0;
                ip = 1;
                j0 = s;
                jp = S; 
            } else {
                i0 = s;
                ip = S;
                j0 = 0;
                jp = 1; 
            }
		
	    
            for(int i=i0; i<I; i+=ip) 
            {
                for(int j=j0; j<J; j+=jp) 
                {		    
                    int x0 = (idw[j]-nb>0) ? (idw[j]-nb) : 0;
                    int y0 = (idh[i]-mb>0) ? (idh[i]-mb) : 0; 
                    int x1 = (nb-idw[j]>0) ? (nb-idw[j]) : 0; 
                    int y1 = (mb-idh[i]>0) ? (mb-idh[i]) : 0; 
                    int M = (idh[i]<ma) ? idh[i]-y0 : ma-y0;
                    int N = (idw[j]<na) ? idw[j]-x0 : na-x0;
		    		  
		    double ff = 0.0;
		    double mf = 0.0;
                    double fm = 0.0;
                    double mmij = 0.0;
		    double w = 0.0;
                    double tmp = 0.0;
		    
                    if (alpha==0.0) {
                        //double c = 1.0/(double)((ma+mb-1)*(na+nb-1));
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {
				w = m1[y0+m][x0+n]*m2[y1+m][x1+n];
                                ff += w*Math.log(Math.abs(f1[y0+m][x0+n] - f2[y1+m][x1+n]));
                                mmij += w;			
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = Math.exp(ff/mmij); 	
                    } else if(alpha==1.0) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {	
				w = m1[y0+m][x0+n]*m2[y1+m][x1+n];
                                ff += w*Math.abs(f1[y0+m][x0+n] - f2[y1+m][x1+n]);				
                                mf += w*Math.abs(f2[y1+m][x1+n]);
                                fm += w*Math.abs(f1[y0+m][x0+n]);
                                mmij += w;   		
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = ff/mmij;
                    } else if(alpha==2.0) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {
				w = m1[y0+m][x0+n]*m2[y1+m][x1+n];
                                ff += w*f1[y0+m][x0+n]*f2[y1+m][x1+n];				
                                mf += w*f2[y1+m][x1+n]*f2[y1+m][x1+n];
                                fm += w*f1[y0+m][x0+n]*f1[y0+m][x0+n];			
                                mmij += w;
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = Math.sqrt((fm-2*ff+mf)/mmij);
                    } else if(alpha==Double.POSITIVE_INFINITY) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {	
				w = (m1[y0+m][x0+n]*m2[y1+m][x1+n]>0)?1.0:0.0;
				mmij += w;				
                                tmp = w*Math.abs(f1[y0+m][x0+n] - f2[y1+m][x1+n]);
                                ff = (tmp>ff)?tmp:ff; 		                                
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = ff;
                    } else if(Double.isNaN(alpha)) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {	
				w =  m1[y0+m][x0+n]*m2[y1+m][x1+n];
                                ff   += w*f1[y0+m][x0+n]*f2[y1+m][x1+n];                                							
                                mmij += w;
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = ff/mmij;
                    } else {				       		    
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {	
				w = m1[y0+m][x0+n]*m2[y1+m][x1+n];	
                                ff += w*mypow(Math.abs(f1[y0+m][x0+n] - f2[y1+m][x1+n]),alpha);
                                mmij += w;			
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = mypow(ff/mmij,1/alpha);			
                    }		    		    
                    mm[idh[i]-1][idw[j]-1] =  mmij / ((ma>mb)?ma:mb) / ((na>nb)?na:nb);;
		    
		}
            }	    
	    return new Integer(s);
        }
    }

    private class Worker2 implements Callable<Integer> 
    {
        private int s;
        private int S;
	
        public Worker2(int s, int S) 
        { 
    	    this.s = s;
            this.S = S;
    	}

        public Integer call() throws Exception
        {
            int ma = f1.length;
            int na = f1[0].length;
            int mb = f1.length;
            int nb = f1[0].length;

    	    int I = idh.length;
            int J = idw.length;
            int i0,ip,j0,jp;

            if(I<J) {
                i0 = 0;
                ip = 1;
                j0 = s;
                jp = S; 
            } else {
                i0 = s;
                ip = S;
                j0 = 0;
                jp = 1; 
            }
		
	    
            for(int i=i0; i<I; i+=ip) 
            {
                for(int j=j0; j<J; j+=jp) 
                {		    
                    int x0 = (idw[j]-nb>0) ? (idw[j]-nb) : 0;
                    int y0 = (idh[i]-mb>0) ? (idh[i]-mb) : 0; 
                    int x1 = (nb-idw[j]>0) ? (nb-idw[j]) : 0; 
                    int y1 = (mb-idh[i]>0) ? (mb-idh[i]) : 0; 
                    int M = (idh[i]<ma) ? idh[i]-y0 : ma-y0;
                    int N = (idw[j]<na) ? idw[j]-x0 : na-x0;
		    		  
		    double ff = 0.0;
		    double mf = 0.0;
                    double fm = 0.0;
                    double mmij = 0.0;
                    double ffm = 0.0;
                    double mff = 0.0;
                    double tmp = 0.0;
		    double w = 0.0;
                    double num = 0.0;
                    double den1 = 0.0;
                    double den2 = 0.0;	
		    
                    if (alpha==0.0) {
                        //double c = 1.0/(double)((ma+mb-1)*(na+nb-1));
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {	
				w = m1[y0+m][x0+n]*m2[y1+m][x1+n];	
                                ff += w*Math.log(Math.abs(f1[y0+m][x0+n] - f2[y1+m][x1+n]));				
                                mf += w*Math.log(Math.abs(f2[y1+m][x1+n]));
                                fm += w*Math.log(Math.abs(f1[y0+m][x0+n]));
                                mmij += w;		
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = Math.exp((ff-.5*mf-.5*fm)/mmij); 	
                    } else if(alpha==1.0) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
				w = m1[y0+m][x0+n]*m2[y1+m][x1+n];
                                ff += w*Math.abs(f1[y0+m][x0+n] - f2[y1+m][x1+n]);				
                                mf += w*Math.abs(f2[y1+m][x1+n]);
                                fm += w*Math.abs(f1[y0+m][x0+n]);
                                mmij += w;			
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = ff/Math.sqrt(fm*mf);	
                    } else if(alpha==2.0) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {
				w = m1[y0+m][x0+n]*m2[y1+m][x1+n];
                                ff += w*f1[y0+m][x0+n]*f2[y1+m][x1+n];				
                                mf += w*f2[y1+m][x1+n]*f2[y1+m][x1+n];
                                fm += w*f1[y0+m][x0+n]*f1[y0+m][x0+n];			
                                mmij += w;
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = Math.sqrt((fm-2*ff+mf)/Math.sqrt(fm*mf));
                    } else if(alpha==Double.POSITIVE_INFINITY) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {	
				w = (m1[y0+m][x0+n]*m2[y1+m][x1+n])>0?1.0:0.0;
                                tmp = w*Math.abs(f1[y0+m][x0+n] - f2[y1+m][x1+n]);
                                ff = (tmp>ff)?tmp:ff; 				
                                tmp = w*Math.abs(f2[y1+m][x1+n]);
                                mf = (tmp>mf)?tmp:mf; 
                                tmp = w*Math.abs(f1[y0+m][x0+n]);	
                                fm = (tmp>fm)?tmp:fm;		
                                mmij += w;
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = ff/Math.sqrt(fm*mf);
                    } else if(Double.isNaN(alpha)) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {	
				w = m1[y0+m][x0+n]*m2[y1+m][x1+n];
                                ff  += w*f1[y0+m][x0+n]*f2[y1+m][x1+n];
                                tmp  = w*f1[y0+m][x0+n];
                                fm  += tmp;			
                                ffm += f1[y0+m][x0+n]*tmp;	
                                tmp  = w*f2[y1+m][x1+n];
                                mf  += tmp;
                                mff += tmp*f2[y1+m][x1+n];							
                                mmij  += w; 
                            }
                        } 
                        num = ff - fm*mf/mmij;
                        den1 = ffm - fm*fm/mmij;
                        den2 = mff - mf*mf/mmij;
                        cc[idh[i]-1][idw[j]-1] = num/Math.sqrt(den1*den2);
                    } else {				       		    
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {	
				w = m1[y0+m][x0+n]*m2[y1+m][x1+n];	
                                ff += w*mypow(Math.abs(f1[y0+m][x0+n] - f2[y1+m][x1+n]),alpha);				
                                mf += w*mypow(Math.abs(f2[y1+m][x1+n]),alpha);
                                fm += w*mypow(Math.abs(f1[y0+m][x0+n]),alpha);
                                mmij += w;		
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = mypow(ff/Math.sqrt(fm*mf),1/alpha);			
                    }		    		    
                    mm[idh[i]-1][idw[j]-1] =  mmij / ((ma>mb)?ma:mb) / ((na>nb)?na:nb);;
		    
		}
            }
	    return new Integer(s);
        }
    }

    public double[][] getmm() 
    {
        return mm;
    }
    
    public double[][] getcc() 
    {
        return cc;
    }

    private double mypow(double x,double alpha) 
    {
        return Math.exp(alpha*Math.log(x));
    }

}

