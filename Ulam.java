import java.util.*;
import java.io.PrintWriter;
import java.io.IOException;



public class Ulam {
    static int maxn = 100000; // Adjust this value for testing purposes
    static int[] a = new int[maxn + 1]; 
    static int[] nx = new int[maxn + 1];
    static int[] pv = new int[maxn + 1]; 
    static int[] k = new int[maxn / 2]; 
    static int nindex = maxn / 100;
    static int[] index = new int[nindex];
    static int nbin = 12000; 
    static int[] bins = new int[nbin]; 
    static long kk1 = 0;
    static long kk2 = 0;
    static long kk3 = 0;
    static long kk4 = 0;
    static long kk5 = 0;
    static double lamda =  2.44344296778474

            ;

    static double step = 13.517831473;
    static List<Long> ulamNumbers = new ArrayList<>();

    public static void main(String[] args) {
        double alpha = 2.0 * Math.PI / lamda;
        double lamdarun = lamda;
        System.out.println("lamda = " + lamda);
        System.out.println("alpha = " + alpha);
        initUlam();
        // initialise index
        for (int i = 0; i < nindex; i++) {
            index[i] = 0;
        }
        pv[0] = 0; 
        nx[0] = 0; 
        setUlam(0, 0);
        setUlam(1, 1);
        setlinks(1);
        setUlam(2, 2);
        setlinks(2);
        int n = 2;
        int nol = 1;
        int nor = 0;
        long bestgap = 0;
        for (long a0 = 3; n < maxn; a0++) {
          
            double rd0 = mod(a0, lamda) / lamda;
            boolean more = true;
            int kount2 = 0;
            long a1x = 0;
            boolean ulam = false;
            if (rd0 < 0.24 || rd0 > 0.80) { 
                int j = n; 
                long aj = getUlam(j);
                while (more && 2 * aj > a0) {
                    kount2++;
                    long a1 = aj;
                    long a2 = a0 - a1;
                    kk3++;
                    if (isUlam(a2)) {
                        if (ulam) { 
                            ulam = false;
                            more = false;
                        } else {
                            ulam = true;
                            a1x = a2;
                        }
                    }
                    j--;
                    aj = getUlam(j);
                }
                more = false;
            }
            int kount0 = 0;
            int i = nx[0]; 
            long ai = getUlam(i);
            double rdi = mod(ai, lamda) / lamda;
            while (2 * rdi <= rd0 + 0.0002 && more && i != 0) {
                kount0++;
                long a2 = a0 - ai;
                kk1++;
                if (isUlam(a2) && ai != a2 && a2 != a1x) { 
                    if (ulam) { 
                        more = false; 
                        ulam = false;
                    } else { 
                        ulam = true;
                        a1x = ai; 
                    }
                }
                i = nx[i]; 
                ai = getUlam(i);
                rdi = mod(ai, lamda) / lamda;
            }
            int kount1 = 0;
            i = pv[0]; 
            ai = getUlam(i);
            rdi = mod(ai, lamda) / lamda;
            while (2 * (1.0 - rdi) <= (1.0 - rd0) + 0.0002 && more && i != 0) {
                kount1++;
                long a2 = a0 - ai;
                kk2++;
                if (isUlam(a2) && ai != a2 && a2 != a1x) { // pair adds up
                    if (ulam) { 
                        more = false;
                        ulam = false;
                    } else { 
                        ulam = true;
                        a1x = ai; 
                    }
                }
                i = pv[i]; 
                ai = getUlam(i);
                rdi = mod(ai, lamda) / lamda;
            }
            if (ulam) {
                n++;
                setUlam(a0, n);
                ulamNumbers.add(a0);
                double z = mod(a0, lamda) / lamda;
                setlinks(n);
                long d = (long) (a0 / lamda);
                double p = 0.0;
                if (z < 1.0 / 3.0) {
                    nor++;
                    p = a0 / (d + 1.0 / 3.0);
                }
                if (z > 2.0 / 3.0) {
                    nol++;
                    p = a0 / (d + 2.0 / 3.0);
                }
                if (z > 2.0 / 3.0 || z < 1.0 / 3.0) {
                    lamdarun = (lamdarun * 9.0 + p) / 10.0;
                    // System.out.println(nor+" "+nol+" "+n+" "+a0+" "+z+" "+p+" "+lamdarun);
                   
                          //  (kk3 / a0) + " " + (kk4 / a0) + " " + (kk5 / a0));
                }
                if (n % 1000000 == 0) {
                    System.out.println("a[" + n + "] = " + a0 + " = " + a1x + " + " + (a0 - a1x));
                    //System.err.println("a[" + n + "] = " + a0 + " = " + a1x + " + " + (a0 - a1x));
                }
                long a1 = getUlam(n - 1);
                long gap = a0 - a1;
                if (gap > bestgap) {
                    bestgap = gap;
                    //System.out.println(n + " " + a0 + " - " + a1 + " = " + gap + " is bigger gap");
                }
                // build distribution by residue in bins
                int ibin = (int) (z * nbin);
                bins[ibin]++;
            }
        }
       
       // System.out.println("biggest gap was " + bestgap);
        double density = ((double) n) / ((double) getUlam(n));
        System.out.println("density = " + density);
        System.out.println("step = " + (1.0 / density));
        System.out.println("");
      
        for (int ibin = 0; ibin < nbin; ibin++) {
            System.out.println(ibin + "," + bins[ibin]);
        }
        checklinks(n);

        // prints  the Sequence to a text file
        try {
            PrintWriter writer = new PrintWriter("Ulam_output.txt");
            for (Long ulamNumber : ulamNumbers) {
                writer.println(ulamNumber);  // This ensures each number is on a new line
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


        // Print Ulam numbers
        System.out.println("First " + maxn + " Ulam numbers:");
        //System.out.println(ulamNumbers);
    }



    public static double mod(long x, double m) {
        double dx = x;
        double z = dx / m;
        long iz = (long) z;
        z -= iz;
        z *= m;
        return z;
    }

    public static void setlinks(int n) {
        double rdn = mod(getUlam(n), lamda) / lamda;
        int j = (int) (nindex * rdn);
        int pvi = index[j];
        boolean more = true;
        while (more) {
            kk4++;
            int i = nx[pvi];
            long ai = getUlam(i);
            double rdi = mod(ai, lamda) / lamda;
            if (i == 0) {
                more = false;
            } else if (rdi < rdn) {
                pvi = i;
            } else {
                more = false;
            }
        }
        int nxi = nx[pvi];
        pv[n] = pvi;
        nx[pvi] = n;
        nx[n] = nxi;
        pv[nxi] = n;
        // update index
        j++;
        double rdi = 0;
        if (j < nindex) rdi = mod(getUlam(index[j]), lamda) / lamda;
        while (j < nindex && rdi < rdn) {
            kk5++;
            index[j] = n;
            j++;
            if (j < nindex) rdi = mod(getUlam(index[j]), lamda) / lamda;
        }
    }

    static void checklinks(int n) {
        int pvi = 0;
        int m = 0;
        while (nx[pvi] != 0) {
            if (pv[nx[pvi]] != pvi)
                System.err.println("links are inconsistent at " + pvi + " -> " + nx[pvi] + " <- " + pv[nx[pvi]]);
            if (nx[pv[pvi]] != pvi)
                System.err.println("links are inconsistent at " + pvi + " <- " + pv[pvi] + " -> " + nx[pv[pvi]]);
            double rdi = mod(getUlam(pvi), lamda) / lamda;
            double rdn = mod(getUlam(nx[pvi]), lamda) / lamda;
            if (rdi > rdn)
                System.err.println("links are not ordered at " + pvi + "," + nx[pvi] + " " + rdi + " > " + rdn);
            m++;
            pvi = nx[pvi];
        }
        if (m != n) System.err.println("link list is wrong length " + m + " != " + n);
        System.err.println("links check complete");
    }

    // booleans flagging the ulam numbers are packed to save space
    static int[] pow2 = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768,
            1 << 16, 1 << 17, 1 << 18, 1 << 19, 1 << 20, 1 << 21, 1 << 22,
            1 << 23, 1 << 24, 1 << 25, 1 << 26, 1 << 27, 1 << 28, 1 << 29};

    public static boolean isUlam(long a0) {
        long i30 = 30;
        int m30 = (int) (a0 % i30);
        int d30 = (int) (a0 / i30);

        if (m30 < 0) {
           m30 += 30; 
        }

        boolean ulam = ((k[d30] & pow2[m30]) > 0);
        return ulam;
    }

    public static void setUlam(long a0, int n) {
        long i30 = 30;
        int m30 = (int) (a0 % i30);
        int d30 = (int) (a0 / i30);
        k[d30] |= pow2[m30];
        double dn = (double) n;
        long ground = (long) (dn * step);
        int an = (int) (a0 - ground);
        a[n] = an;
    }

    public static long getUlam(int n) {
        double dn = (double) n;
        long ground = (long) (dn * step);
        long a0 = ground + a[n];
        return a0;
    }

    public static void initUlam() {
        Arrays.fill(k, 0);
        Arrays.fill(bins, 0);
    }
}
