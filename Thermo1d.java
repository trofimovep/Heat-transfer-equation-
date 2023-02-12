public class Thermo1d {

    public static void main(String[] args) {
        System.out.println("Hello!");
    }
    
    public static final double[] heat1D(
            int cond, double rho, double c, double t, double length, double h,
            double T0, double lambda, double alpha, double Tleft, double Tright
    ) {
        
        double a = lambda / (rho * c);
        double tau = (h * h) / (4 * a);
        int Nt = (int) Math.round(t / tau);
        int Nl = (int) Math.round(length / h + 1);
        
        double[][] T = new double[Nl][Nt];
        for (int i = 0; i < Nl; i++) {
            T[i][0] = T0;
        }
        if (cond == 1) {
            for (int j = 0; j < Nt - 1; j++) {
                for (int i = 0; i < Nl; i++) {
                    if (i == 0) {
                        T[i][j + 1] = Tleft;
                    } else if (i == Nl - 1) {
                        T[i][j + 1] = Tright;
                    } else
                        T[i][j + 1] = T[i][j] + a * tau * (T[i - 1][j] - 2 * T[i][j] + T[i + 1][j]) / (h * h);
                }
            }
        }
        if (cond == 2) {
            for (int j = 0; j < Nt - 1; j++) {
                for (int i = 0; i < Nl; i++) {
                    if (i == 0) {
                        T[i][j + 1] = T[i][j] + (Tleft + ((lambda * (T[i + 1][j] - T[i][j])) / h)) * (2 * tau / (c * rho * h));
                    } else if (i == Nl - 1) {
                        T[i][j + 1] = T[i][j] + (Tleft + ((lambda * (T[i - 1][j] - T[i][j])) / h)) * (2 * tau / (c * rho * h));
                    } else {
                        T[i][j + 1] = T[i][j] + tau * a * (T[i - 1][j] - 2 * T[i][j] + T[i + 1][j]) / (h * h);
                    }
                }
            }
        }
        if (cond == 3) {
            for (int j = 0; j < Nt - 1; j++) {
                for (int i = 0; i < Nl; i++) {
                    if (i == 0) {
                        T[i][j + 1] = T[i][j] + (((alpha * (Tleft - T[i][j])) + ((lambda * (T[i + 1][j] - T[i][j]))) / h)) * (2 * tau / (c * rho * h));
                    } else if (i == Nl - 1) {
                        T[i][j + 1] = T[i][j] + (((alpha * (Tright - T[i][j])) + ((lambda * (T[i - 1][j] - T[i][j]))) / h)) * (2 * tau / (c * rho * h));
                        
                    } else
                        T[i][j + 1] = T[i][j] + tau * a * (T[i - 1][j] - 2 * T[i][j] + T[i + 1][j]) / (h * h);
                }
            }
        }
        double[] Tend = new double[Nl];
        for (int i = 0; i < Nl; i++)
            Tend[i] = T[i][Nt - 1];
        for (int i = 0; i < Nl; i++)
            System.out.println("Tend = " + Tend[i]);
        System.out.println("alpha = " + alpha);
        return Tend;
    }
}
