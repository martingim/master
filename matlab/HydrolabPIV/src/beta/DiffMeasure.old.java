import java.lang.Thread;

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
    
    private static int maxthreads = Integer.MAX_VALUE;

    public DiffMeasure() 
    {
       
    }
    
    public DiffMeasure(int n) 
    {
        maxthreads = n;
    }

    public double[][] maskedncc(double[][] f1, double[][] m1, double[][] f2, double[][] m2, int[] idh, int[] idw) 
    {
        return maskednmd(f1,m1,f2,m2,idh,idw,Double.NaN);
    }

    public double[][] maskedncc(double[][] f1, double[][] m1, double[][] f2, double[][] m2, int[] idh, int[] idw, double alpha) 
    {        
        return maskednmd(f1,m1,f2,m2,idh,idw,Double.NaN);
    }

    public double[][] maskednmd(double[][] f1, double[][] m1, double[][] f2, double[][] m2, int[] idh, int[] idw, double alpha) 
    {	
        int ma = f1.length;
        int na = f1[0].length;
        int mb = f1.length;
        int nb = f1[0].length;
        
        for(int m=0; m<ma; m++) {
    	    for(int n=0; n<na; n++) {
                f1[m][n] *= m1[m][n];
            }
        }
        for(int m=0; m<mb; m++) {
            for(int n=0; n<nb; n++) {
                f2[m][n] *= m2[m][n];
            }
        }

        this.cc = new double[ma+mb-1][na+nb-1];
        this.mm = new double[ma+mb-1][na+nb-1];
		
    	int S =  Runtime.getRuntime().availableProcessors(); // number of threads
    	S = (S<=maxthreads) ? S : maxthreads;
    
    	if(S==1) {
            for(int i=0; i<idh.length; i++) {
                for(int j=0; j<idw.length; j++) {    
                    int x0 = (idw[j]-nb>0) ? (idw[j]-nb) : 0;
                    int y0 = (idh[i]-mb>0) ? (idh[i]-mb) : 0; 
                    int x1 = (nb-idw[j]>0) ? (nb-idw[j]) : 0; 
                    int y1 = (mb-idh[i]>0) ? (mb-idh[i]) : 0; 
                    int M = (idh[i]<ma) ? idh[i]-y0 : ma-y0;
                    int N = (idw[j]<na) ? idw[j]-x0 : na-x0;
		    
        		    double ff = 0.0;
                    double mf = 0.0;
                    double fm = 0.0;
                    double mm = 0.0;
                    double ffm = 0.0;
                    double mff = 0.0;
                    double tmp = 0.0;
                    double num = 0.0;
                    double den1 = 0.0;
                    double den2 = 0.0;		   
    
        		    if (alpha==0.0) {
                    	//double c = 1.0/(double)((ma+mb-1)*(na+nb-1));
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
                                ff += Math.log(Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n] - m1[y0+m][x0+n]*f2[y1+m][x1+n]));				
                                mf += Math.log(Math.abs(m1[y0+m][x0+n]*f2[y1+m][x1+n]));
                                fm += Math.log(Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n]));
                                mm += m1[y0+m][x0+n]*m2[y1+m][x1+n];			
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = Math.exp((ff-.5*mf-.5*fm)/mm); 	
                    } else if(alpha==1.0) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
                                ff += Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n] - m1[y0+m][x0+n]*f2[y1+m][x1+n]);				
                                mf += Math.abs(m1[y0+m][x0+n]*f2[y1+m][x1+n]);
                                fm += Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n]);
                                mm += m1[y0+m][x0+n]*m2[y1+m][x1+n];			
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = ff/Math.sqrt(fm*mf);	
                    } else if(alpha==2.0) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
                                ff += f1[y0+m][x0+n]*f2[y1+m][x1+n];				
                                mf += m1[y0+m][x0+n]*f2[y1+m][x1+n]*f2[y1+m][x1+n];
                                fm += f1[y0+m][x0+n]*f1[y0+m][x0+n]*m2[y1+m][x1+n];			
                                mm += m1[y0+m][x0+n]*m2[y1+m][x1+n];
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = Math.sqrt((fm-2*ff+mf)/Math.sqrt(fm*mf));
                    } else if(alpha==Double.POSITIVE_INFINITY) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
                                tmp = Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n] - m1[y0+m][x0+n]*f2[y1+m][x1+n]);
                                ff = (tmp>ff)?tmp:ff; 				
                                tmp = Math.abs(m1[y0+m][x0+n]*f2[y1+m][x1+n]);
                                mf = (tmp>mf)?tmp:mf; 
                                tmp = Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n]);	
                                fm = (tmp>fm)?tmp:fm;		
                                mm += m1[y0+m][x0+n]*m2[y1+m][x1+n];
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = ff/Math.sqrt(fm*mf);
                    } else if(Double.isNaN(alpha)) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
                                ff  += f1[y0+m][x0+n]*f2[y1+m][x1+n];
                                tmp  = f1[y0+m][x0+n]*m2[y1+m][x1+n];
                                fm  += tmp;			
                                ffm += f1[y0+m][x0+n]*tmp;	
                                tmp  = m1[y0+m][x0+n]*f2[y1+m][x1+n];
                                mf  += tmp;
                                mff += tmp*f2[y1+m][x1+n];							
                                mm  += m1[y0+m][x0+n]*m2[y1+m][x1+n];
                            }
                        } 
                        num = ff - fm*mf/mm;
                        den1 = ffm - fm*fm/mm;
                        den2 = mff - mf*mf/mm;
                        cc[idh[i]-1][idw[j]-1] = num/Math.sqrt(den1*den2);
                    } else {				       		    
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
                                ff += mypow(Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n] - m1[y0+m][x0+n]*f2[y1+m][x1+n]),alpha);				
                                mf +=  mypow(Math.abs(m1[y0+m][x0+n]*f2[y1+m][x1+n]),alpha);
                                fm +=  mypow(Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n]),alpha);
                                mm += m1[y0+m][x0+n]*m2[y1+m][x1+n];			
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = mypow(ff/Math.sqrt(fm*mf),1/alpha);			
                    }
    
		    
            	    this.mm[idh[i]-1][idw[j]-1] =  mm / ((ma>mb)?ma:mb) / ((na>nb)?na:nb);
		
                }
            }
        } else {
            this.f1 = f1;
            this.m1 = m1;
            this.f2 = f2;
            this.m2 = m2;
            this.idh = idh;
            this.idw = idw;
            this.alpha = alpha;
    
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
        return cc;
    }
       
    private class Worker implements Runnable 
    {
        private int s;
        private int S;
	
        public Worker(int s, int S) 
        { 
    	    this.s = s;
            this.S = S;
    	}

        public void run() 
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
                    double num = 0.0;
                    double den1 = 0.0;
                    double den2 = 0.0;	
		    
                    if (alpha==0.0) {
                        //double c = 1.0/(double)((ma+mb-1)*(na+nb-1));
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
                                ff += Math.log(Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n] - m1[y0+m][x0+n]*f2[y1+m][x1+n]));				
                                mf += Math.log(Math.abs(m1[y0+m][x0+n]*f2[y1+m][x1+n]));
                                fm += Math.log(Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n]));
                                mmij += m1[y0+m][x0+n]*m2[y1+m][x1+n];			
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = Math.exp((ff-.5*mf-.5*fm)/mmij); 	
                    } else if(alpha==1.0) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
                                ff += Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n] - m1[y0+m][x0+n]*f2[y1+m][x1+n]);				
                                mf += Math.abs(m1[y0+m][x0+n]*f2[y1+m][x1+n]);
                                fm += Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n]);
                                mmij += m1[y0+m][x0+n]*m2[y1+m][x1+n];			
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = ff/Math.sqrt(fm*mf);	
                    } else if(alpha==2.0) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
                                ff += f1[y0+m][x0+n]*f2[y1+m][x1+n];				
                                mf += m1[y0+m][x0+n]*f2[y1+m][x1+n]*f2[y1+m][x1+n];
                                fm += f1[y0+m][x0+n]*f1[y0+m][x0+n]*m2[y1+m][x1+n];			
                                mmij += m1[y0+m][x0+n]*m2[y1+m][x1+n];
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = Math.sqrt((fm-2*ff+mf)/Math.sqrt(fm*mf));
                    } else if(alpha==Double.POSITIVE_INFINITY) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
                                tmp = Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n] - m1[y0+m][x0+n]*f2[y1+m][x1+n]);
                                ff = (tmp>ff)?tmp:ff; 				
                                tmp = Math.abs(m1[y0+m][x0+n]*f2[y1+m][x1+n]);
                                mf = (tmp>mf)?tmp:mf; 
                                tmp = Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n]);	
                                fm = (tmp>fm)?tmp:fm;		
                                mmij += m1[y0+m][x0+n]*m2[y1+m][x1+n];
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = ff/Math.sqrt(fm*mf);
                    } else if(Double.isNaN(alpha)) {
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
                                ff  += f1[y0+m][x0+n]*f2[y1+m][x1+n];
                                tmp  = f1[y0+m][x0+n]*m2[y1+m][x1+n];
                                fm  += tmp;			
                                ffm += f1[y0+m][x0+n]*tmp;	
                                tmp  = m1[y0+m][x0+n]*f2[y1+m][x1+n];
                                mf  += tmp;
                                mff += tmp*f2[y1+m][x1+n];							
                                mmij  += m1[y0+m][x0+n]*m2[y1+m][x1+n];
                            }
                        } 
                        num = ff - fm*mf/mmij;
                        den1 = ffm - fm*fm/mmij;
                        den2 = mff - mf*mf/mmij;
                        cc[idh[i]-1][idw[j]-1] = num/Math.sqrt(den1*den2);
                    } else {				       		    
                        for(int m=0; m<M; m++) {
                            for(int n=0; n<N; n++) {			
                                ff += mypow(Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n] - m1[y0+m][x0+n]*f2[y1+m][x1+n]),alpha);				
                                mf +=  mypow(Math.abs(m1[y0+m][x0+n]*f2[y1+m][x1+n]),alpha);
                                fm +=  mypow(Math.abs(f1[y0+m][x0+n]*m2[y1+m][x1+n]),alpha);
                                mmij += m1[y0+m][x0+n]*m2[y1+m][x1+n];			
                            }
                        } 
                        cc[idh[i]-1][idw[j]-1] = mypow(ff/Math.sqrt(fm*mf),1/alpha);			
                    }		    		    
                    mm[idh[i]-1][idw[j]-1] =  mmij / ((ma>mb)?ma:mb) / ((na>nb)?na:nb);;
		    
        		}
            }
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

