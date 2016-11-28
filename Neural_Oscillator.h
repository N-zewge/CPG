/* parameter
	tc[0], tc[1] : tau, tau'
	ue, uf, ve, vf : u^e, u^f, v^e, v^f
	ye, yf : y^e, y^f
	wee, wef, wfe, wff : connecting weight
	uinput : u_0 (input signal)
	Neural_N : number of neural oscillators
*/

#include<algorithm>

#define Neural_N 3

double * Neural_Oscillator(double *uee, double *uff, double *vee, double *vff, double *yee, double *yff, double *weee, double *weff, double *wfee, double *wfff, double uinput1, double *tc1)
{
	double Beta = 2.5;
	int i, j;
	static double output_value[Neural_N];
	double due[Neural_N], dve[Neural_N], duf[Neural_N], dvf[Neural_N]; //derivative of u^e, v^e, u^f, v^f
	double Sume[Neural_N], Sumf[Neural_N];
	double Ts = 0.005; //sampling time

	for (i = 0; i < Neural_N; i++)
	{
		Sume[i] = -2.5 * yff[i];
		Sumf[i] = -2.5 * yee[i];
		for (j = 0; j < Neural_N; j++)
		{
			if (j != i)
			{
				Sume[i] -= weee[i * Neural_N + j] * yee[j] + weff[i * Neural_N + j] * yff[j];
				Sumf[i] -= wfee[i * Neural_N + j] * yee[j] + wfff[i * Neural_N + j] * yff[j];
			}
		}
	}
		for (i = 0; i < Neural_N; i++)
		{
			due[i] = (-uee[i] + Sume[i] - Beta * vee[i] + uinput1) / tc1[0];
			dve[i] = (-vee[i] + yee[i]) / tc1[1];
			duf[i] = (-uff[i] + Sumf[i] - Beta * vff[i] + uinput1) / tc1[0];
			dvf[i] = (-vff[i] + yff[i]) / tc1[1];
		}


		for (i = 0; i < Neural_N; i++)
		{
			uee[i] = uee[i] + due[i] * Ts; // current ue is (previous ue + derivative of ue*sampling time)
			vee[i] = vee[i] + dve[i] * Ts;
			yee[i] = fmax(uee[i], 0);

			uff[i] = uff[i] + duf[i] * Ts;
			vff[i] = vff[i] + dvf[i] * Ts;
			yff[i] = fmax(uff[i], 0);
			
			output_value[i] = yee[i] - yff[i];
			
		}
		return output_value;
	}
	
