package mono;

import java.util.Arrays;

public class Mono {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double param1 = 0.15;
		double param2 = 0.1;
		double param3 = 0.6;
		double param4 = 0.3;
		double param5 = 0.5;
		double param6 = 6;
		double param7 = 0.84;
		double param8 = 20;
		double param9 = 1030;	
		double param10 = 0.01;
		double param11 = 1;
		double[] input = {param1,param2,param3,param4,param5,param6,param7,param8,param9,param10,param11};
		double output = calFilter(input);
		System.out.print(output);
	}
	private static double calFilter(double[] input) {
		double col_dia = input[0];
		double in_conc = input[1];
		double media_size = input[2]*1e-3;
		double media_depth = input[3];
		double media_por = input[4];
		double in_vS = input[5];
		double in_mu = input[6]*1e-6;
		double in_nu = in_mu*1e3;
		double in_dia_p = input[7]*1e-6;
		double in_rho_p = input[8];
		double col_deltaz = input[9];
		int fil_time = (int) Math.round(input[10]*3600);
		double ratio = media_size/in_dia_p;
		double SA_filter = (col_dia*col_dia)*Math.PI/4;
		int no_depths_an = (int) Math.round(media_depth/col_deltaz);
		double[][] an_c = new double[no_depths_an][fil_time];
		double[][] an_dc = new double[no_depths_an][fil_time];
		double[][] an_poro = new double[no_depths_an][fil_time];
		for (double[] row: an_poro)
		    Arrays.fill(row, media_por);

		double[][] an_eta_interception = new double[no_depths_an][fil_time];
		double[][] an_eta_sedimentation = new double[no_depths_an][fil_time];
		double[][] an_eta = new double[no_depths_an][fil_time];
		double[][] an_lamda = new double[no_depths_an][fil_time];
		double[][] an_size_var = new double[no_depths_an][fil_time];
		for (double[] row: an_size_var)
		    Arrays.fill(row, media_size);

		double[][] an_sigma = new double[no_depths_an][fil_time];
		double[][] an_hl_1 = new double[no_depths_an][fil_time];
		double[][] an_hl_2 = new double[no_depths_an][fil_time];
		double[][] headloss = new double[no_depths_an][fil_time];
		for (int t=0;t<fil_time;t++) {
			if (t==0) {
				for (int z=0;z<no_depths_an;z++) {
					an_c[z][t] = in_conc;
					an_size_var[z][t]=media_size;
					an_poro[z][t] = media_por;
					an_hl_1[z][t] = 150*in_mu*(in_vS/3600)*Math.pow(1-an_poro[z][t], 2)*(z+1)*col_deltaz/(9.81*Math.pow(an_size_var[z][t], 2)*Math.pow(an_poro[z][t], 3));
					an_hl_2[z][t] = 1.75*Math.pow(in_vS/3600,2)*(1-an_poro[z][t])*(z+1)*col_deltaz/(9.81*an_size_var[z][t]*Math.pow(an_poro[z][t], 3));
					headloss[z][t] = an_hl_1[z][t] + an_hl_2[z][t];
				}
			}else {
				for (int z=0;z<no_depths_an;z++) {
					an_dc[z][t] = (in_conc-an_c[z][t-1])*(t+1)*(in_vS/3600)/(in_rho_p*(z+1)*col_deltaz*media_por);
					an_sigma[z][t] = (in_conc-an_c[z][t-1])*(t+1)*(in_vS/3600)*SA_filter;
					an_poro[z][t]  = media_por - an_dc[z][t];
					an_size_var[z][t] = media_size*(1+an_dc[z][t]);
					an_eta_interception[z][t] = 1.5*Math.pow(in_dia_p/an_size_var[z][t],2) ;
					an_eta_sedimentation[z][t]  = Math.pow(in_dia_p, 2)*9.81*(in_rho_p - 1000)/(18*in_nu*in_vS);
					an_eta[z][t]  = an_eta_interception[z][t] + an_eta_sedimentation[z][t];
					an_lamda[z][t]  = 1.5*(1-an_poro[z][t])*an_eta[z][t] /an_size_var[z][t] ;
					an_c[z][t] = in_conc*Math.exp(-an_lamda[z][t]*(z+1)*col_deltaz);
					an_hl_1[z][t] = 150*in_mu*(in_vS/3600)*Math.pow(1-an_poro[z][t], 2)*(z+1)*col_deltaz/(9.81*Math.pow(an_size_var[z][t], 2)*Math.pow(an_poro[z][t], 3));
					an_hl_2[z][t] = 1.75*Math.pow(in_vS/3600,2)*(1-an_poro[z][t])*(z+1)*col_deltaz/(9.81*an_size_var[z][t]*Math.pow(an_poro[z][t], 3));
					headloss[z][t] = (an_hl_1[z][t] + an_hl_2[z][t] - headloss[z][0])*ratio/((z+1)*col_deltaz);
				}
			}
		}
		double output = headloss[no_depths_an-1][fil_time-1];
		return output;
	}
}
