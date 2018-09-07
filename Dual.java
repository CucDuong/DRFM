package mono;

import java.util.Arrays;

public class Dual {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double param1 = 0.15;
		double param2 = 0.1;
		double param3 = 0.6;
		double param4 = 0.3;
		double param5 = 0.5;
		double param6 = 0.2;
		double param7 = 0.3;
		double param8 = 0.4;
		double param9 = 6.0;	
		double param10 = 0.84;
		double param11 = 20;
		double param12 = 1030;
		double param13 = 0.01;
		double param14 = 1.0;
		double[] input = {param1,param2,param3,param4,param5,param6,param7,param8,param9,param10,param11,param12,param13,param14};
		double[] output = calFilter(input);
		System.out.print(output[0] + " and " + output[1]);
	}
	private static double[] calFilter(double[] input) {
		double col_dia = input[0];
		double in_conc = input[1];
		double media_size1 = input[2]*1e-3;
		double media_depth1 = input[3];
		double media_por1 = input[4];
		double media_size2 = input[5]*1e-3;
		double media_depth2 = input[6];
		double media_por2 = input[7];
		double in_vS = input[8];
		double in_mu = input[9]*1e-6;
		double in_nu = in_mu*1e3;
		double in_dia_p = input[10]*1e-6;
		double in_rho_p = input[11];
		double col_deltaz = input[12];
		int fil_time = (int) Math.round(input[13]*3600);
		double SA_filter = (col_dia*col_dia)*Math.PI/4;
		double ratio1 = media_size1/in_dia_p;
		double ratio2 = media_size2/in_dia_p;
		int no_depths_an = (int) Math.round(media_depth1/col_deltaz);
		int no_depths_sand = (int) Math.round(media_depth2/col_deltaz);
		double[][] an_c = new double[no_depths_an][fil_time];
		double[][] an_dc = new double[no_depths_an][fil_time];
		double[][] an_poro = new double[no_depths_an][fil_time];
		for (double[] row: an_poro)
		    Arrays.fill(row, media_por1);
		double[][] an_eta_interception = new double[no_depths_an][fil_time];
		double[][] an_eta_sedimentation = new double[no_depths_an][fil_time];
		double[][] an_eta = new double[no_depths_an][fil_time];
		double[][] an_lamda = new double[no_depths_an][fil_time];
		double[][] an_size_var = new double[no_depths_an][fil_time];
		for (double[] row: an_size_var)
		    Arrays.fill(row, media_size1);
		double[][] an_sigma = new double[no_depths_an][fil_time];
		double[][] an_hl_1 = new double[no_depths_an][fil_time];
		double[][] an_hl_2 = new double[no_depths_an][fil_time];
		double[][] an_hl_total_metres = new double[no_depths_an][fil_time];
		double[][] headloss_an = new double[no_depths_an][fil_time];
		
		double[][] sand_c= new double[no_depths_sand][fil_time];
		double[][] sand_dc = new double[no_depths_sand][fil_time];
		double[][] sand_poro = new double[no_depths_sand][fil_time];
		for (double[] row: sand_poro)
		    Arrays.fill(row, media_por2);
		double[][] sand_eta_interception = new double[no_depths_sand][fil_time];
		double[][] sand_eta_sedimentation = new double[no_depths_sand][fil_time];
		double[][] sand_eta = new double[no_depths_sand][fil_time];
		double[][] sand_lamda = new double[no_depths_sand][fil_time];
		double[][] sand_size_var = new double[no_depths_sand][fil_time];
		for (double[] row: sand_size_var)
		    Arrays.fill(row, media_size2);
		double[][] sand_sigma = new double[no_depths_sand][fil_time];
		double[][] sand_hl_1 = new double[no_depths_sand][fil_time];
		double[][] sand_hl_2 = new double[no_depths_sand][fil_time];
		double[][] sand_hl_total_metres = new double[no_depths_sand][fil_time];
		double[][] headloss_sand = new double[no_depths_sand][fil_time];
		
		double sand_conc = 0;
		
		for (int t=0;t<fil_time;t++) {
			if (t==0) {
				for (int z=0;z<no_depths_an;z++) {
					
					// for the first one
					an_c[z][t] = in_conc;
					an_size_var[z][t]=media_size1;
					an_poro[z][t] = media_por1;
					an_hl_1[z][t] = 150*in_mu*(in_vS/3600)*Math.pow(1-an_poro[z][t], 2)*(z+1)*col_deltaz/(9.81*Math.pow(an_size_var[z][t], 2)*Math.pow(an_poro[z][t], 3));
					an_hl_2[z][t] = 1.75*Math.pow(in_vS/3600,2)*(1-an_poro[z][t])*(z+1)*col_deltaz/(9.81*an_size_var[z][t]*Math.pow(an_poro[z][t], 3));
					an_hl_total_metres[z][t]=an_hl_1[z][t] + an_hl_2[z][t];
					headloss_an[z][t] = (an_hl_1[z][t] + an_hl_2[z][t])/((z+1)*col_deltaz) ;
				}
				
				for (int z=0;z<no_depths_sand;z++) {	
					// for the second one
					sand_c[z][t] = in_conc;
					sand_size_var[z][t]=media_size2;
					sand_poro[z][t] = media_por2;
					sand_hl_1[z][t] = 150*in_mu*(in_vS/3600)*Math.pow(1-sand_poro[z][t], 2)*(z+1)*col_deltaz/(9.81*Math.pow(sand_size_var[z][t], 2)*Math.pow(sand_poro[z][t], 3));
					sand_hl_2[z][t] = 1.75*Math.pow(in_vS/3600,2)*(1-sand_poro[z][t])*(z+1)*col_deltaz/(9.81*sand_size_var[z][t]*Math.pow(sand_poro[z][t], 3));
					sand_hl_total_metres[z][t] = sand_hl_1[z][t] + sand_hl_2[z][t];
					headloss_sand[z][t] = (sand_hl_1[z][t] + sand_hl_2[z][t])/((z+1)*col_deltaz) ;
				}
			}else {
				for (int z=0;z<no_depths_an;z++) {
					
					// for the first one
					an_dc[z][t] = (in_conc-an_c[z][t-1])*(t+1)*(in_vS/3600)/(in_rho_p*(z+1)*col_deltaz*media_por1);
					an_sigma[z][t] = (in_conc-an_c[z][t-1])*(t+1)*(in_vS/3600)*SA_filter;
					an_poro[z][t]  = media_por1 - an_dc[z][t];
					an_size_var[z][t] = media_size1*(1+an_dc[z][t]);
					an_eta_interception[z][t] = 1.5*Math.pow(in_dia_p/an_size_var[z][t],2) ;
					an_eta_sedimentation[z][t]  = Math.pow(in_dia_p, 2)*9.81*(in_rho_p - 1000)/(18*in_nu*in_vS);
					an_eta[z][t]  = an_eta_interception[z][t] + an_eta_sedimentation[z][t];
					an_lamda[z][t]  = 1.5*(1-an_poro[z][t])*an_eta[z][t] /an_size_var[z][t] ;
					an_c[z][t] = in_conc*Math.exp(-an_lamda[z][t]*(z+1)*col_deltaz);
					
					an_hl_1[z][t] = 150*in_mu*(in_vS/3600)*Math.pow(1-an_poro[z][t], 2)*(z+1)*col_deltaz/(9.81*Math.pow(an_size_var[z][t], 2)*Math.pow(an_poro[z][t], 3));
					an_hl_2[z][t] = 1.75*Math.pow(in_vS/3600,2)*(1-an_poro[z][t])*(z+1)*col_deltaz/(9.81*an_size_var[z][t]*Math.pow(an_poro[z][t], 3));
					headloss_an[z][t] = (an_hl_1[z][t] + an_hl_2[z][t] - an_hl_total_metres[z][0])*ratio1/((z+1)*col_deltaz);
				}
				sand_conc = an_c[no_depths_an-1][t];
				
				for (int z=0;z<no_depths_sand;z++) {
	
					// for the second one
					sand_dc[z][t] = (sand_conc-sand_c[z][t-1])*(t+1)*(in_vS/3600)/(in_rho_p*(z+1)*col_deltaz*media_por2);
					sand_sigma[z][t] = (in_conc-sand_c[z][t-1])*(t+1)*(in_vS/3600)*SA_filter;
					sand_poro[z][t]  = media_por2 - sand_dc[z][t];
					sand_size_var[z][t] = media_size2*(1+sand_dc[z][t]);
					sand_eta_interception[z][t] = 1.5*Math.pow(in_dia_p/sand_size_var[z][t],2) ;
					sand_eta_sedimentation[z][t]  = Math.pow(in_dia_p, 2)*9.81*(in_rho_p - 1000)/(18*in_nu*in_vS);
					sand_eta[z][t]  = sand_eta_interception[z][t] + sand_eta_sedimentation[z][t];
					sand_lamda[z][t]  = 1.5*(1-sand_poro[z][t])*sand_eta[z][t] /sand_size_var[z][t] ;
					sand_c[z][t] = sand_conc*Math.exp(-sand_lamda[z][t]*(z+1)*col_deltaz);
					
					sand_hl_1[z][t] = 150*in_mu*(in_vS/3600)*Math.pow(1-sand_poro[z][t], 2)*(z+1)*col_deltaz/(9.81*Math.pow(sand_size_var[z][t], 2)*Math.pow(sand_poro[z][t], 3));
					sand_hl_2[z][t] = 1.75*Math.pow(in_vS/3600,2)*(1-sand_poro[z][t])*(z+1)*col_deltaz/(9.81*sand_size_var[z][t]*Math.pow(sand_poro[z][t], 3));
					headloss_sand[z][t] = (sand_hl_1[z][t] + sand_hl_2[z][t] - sand_hl_total_metres[z][0])*ratio2/((z+1)*col_deltaz);
					
				}
			}
		}
		double[] output = new double[2];
		output[0] = headloss_an[no_depths_an-1][fil_time-1];
		output[1] = headloss_sand[no_depths_sand-2][fil_time-2];
		return output;
	}

}
