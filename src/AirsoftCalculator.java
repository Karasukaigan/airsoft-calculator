/** 
* A calculator for Airsoft.
* 
* 
* @author  Karasukaigan 
* @date    2021/4/16 
* @version 1.0.0
*/

import java.math.BigDecimal;

public class AirsoftCalculator {
	
	/**
	 * Used to calculate muzzle kinetic energy.
	 * 銃口での運動エネルギーを計算するために使用されます。
	 * 
	 * @param bulletWeight
	 *        Bullet weight. The unit is g.弾丸の重さ。単位はgです。
	 * @param velocity
	 *        The speed of the bullet at the muzzle. The unit is m/s.銃口での弾丸の速度。単位はm/sです。
	 * @return Kinetic energy 弾丸の運動エネルギー、単位はJです。
	 */
	public double calculateKineticEnergy(double bulletWeight , double velocity) {
		double e = 0.0 , m = bulletWeight/1000 , v = velocity , kineticEnergy = e;
		e = 0.5*m*v*v; //unit:J=kg*m^2/s^2
		kineticEnergy = e;
		return kineticEnergy;
	}
	
	/**
	 * Used to calculate muzzle kinetic energy.
	 * Keep only two decimal places.
	 * 銃口での運動エネルギーを計算するために使用されます。
	 * 小数点以下3桁を保持します。
	 * 
	 * @param bulletWeight
	 *        Bullet weight. The unit is g.弾丸の重さ。単位はgです。
	 * @param velocity
	 *        The speed of the bullet at the muzzle. The unit is m/s.銃口での弾丸の速度。単位はm/sです。
	 * @return Kinetic energy 弾丸の運動エネルギー、単位はJです。
	 */	
	public double calculateKineticEnergyRounding(double bulletWeight , double velocity) {
		double e = 0.0 , m = bulletWeight/1000 , v = velocity , kineticEnergy = e;
		e = 0.5*m*v*v; //unit:J=kg*m^2/s^2
		kineticEnergy = (double)Math.round(e * 1000) / 1000;
		return kineticEnergy;
	}

	/**
	 * Used to calculate muzzle kinetic energy.
	 * You can set how many decimal places to keep.
	 * 銃口での運動エネルギーを計算するために使用されます。
	 * 保持したい小数点以下の桁数が設定できます。
	 * 
	 * @param bulletWeight
	 *        Bullet weight. The unit is g.弾丸の重さ。単位はgです。
	 * @param velocity
	 *        The speed of the bullet at the muzzle. The unit is m/s.銃口での弾丸の速度。単位はm/sです。
	 * @param significantFigures(1~)
	 *        保持したい小数点以下の桁数です。（1から）
	 * @return Kinetic energy 弾丸の運動エネルギー、単位はJです。
	 */
	@SuppressWarnings("deprecation")
	public double calculateKineticEnergyRounding(double bulletWeight , double velocity , int significantFigures) {
		double e = 0.0 , m = bulletWeight/1000 , v = velocity , kineticEnergy = e;
		e = 0.5*m*v*v; //unit:J=kg*m^2/s^2
		BigDecimal bd = new BigDecimal(e);
		kineticEnergy = bd.setScale(significantFigures,BigDecimal.ROUND_HALF_UP) .doubleValue();
		return kineticEnergy;
	}
	
	/**
	 * Used to calculate specific kinetic energy (Chinese standard).
	 * 用于计算比动能（中国标准J/cm^2）。
	 * 
	 * @param bulletWeight
	 *        子弹重量。单位是g。
	 * @param velocity
	 *        枪口初速。单位是m/s。
	 * @param bulletDiameter
	 *        子弹直径。单位是mm。
	 * @return specific kinetic energy 比动能
	 */
	public double calculateSpecificKineticEnergyChina(double bulletWeight , double velocity,  double bulletDiameter) {
		double e = 0.0 , m = bulletWeight/1000 , v = velocity , kineticEnergy = e;
		e = 0.5*m*v*v; //unit:J=kg*m^2/s^2
		kineticEnergy = e;
		double pi = 3.1415926535898;
		double a = (bulletDiameter/10/2) * (bulletDiameter/10/2) * pi; //unit:cm^2
		double specificKineticEnergy = (double)Math.round((kineticEnergy / a) * 100) / 100; //unit:J/cm^2
		return specificKineticEnergy;
	}
	
	/**
	 * Used to calculate specific kinetic energy (Chinese standard).
	 * The bullet is a unique 7.5mm in diameter.
	 * 用动能换算比动能（中国标准J/cm^2）。
	 * 7.5mm水弹专用。
	 * 
	 * @param kineticEnergy
	 *        枪口动能。单位是J。
	 * @return specific kinetic energy 比动能
	 */
	public double toSpecificKineticEnergyChina(double kineticEnergy) {
		double bulletDiameter = 7.5; //unit:mm
		double pi = 3.1415926535898;
		double a = (bulletDiameter/10/2) * (bulletDiameter/10/2) * pi; //unit:cm^2
		double specificKineticEnergy = (double)Math.round((kineticEnergy / a) * 100) / 100; //unit:J/cm^2
		return specificKineticEnergy;
	}

	/**
	 * Used to calculate specific kinetic energy (Chinese standard).
	 * 用动能换算比动能（中国标准J/cm^2）。
	 * 可自行设定子弹直径。
	 * 
	 * @param kineticEnergy
	 *        枪口动能。单位是J。
	 * @param bulletDiameter
	 *        子弹直径。单位是mm。
	 * @return specific kinetic energy 比动能
	 */
	public double toSpecificKineticEnergyChina(double kineticEnergy , double bulletDiameter) {
		double pi = 3.1415926535898;
		double a = (bulletDiameter/10/2) * (bulletDiameter/10/2) * pi; //unit:cm^2
		double specificKineticEnergy = (double)Math.round((kineticEnergy / a) * 100) / 100; //unit:J/cm^2
		return specificKineticEnergy;
	}
	
	/**
	 * Check whether the kinetic energy is legal. 
	 * You can set the upper limit of kinetic energy.
	 * If "true" is returned, it is legal.
	 * 運動エネルギーが合法かどうかをチェックします。
	 * 運動エネルギーの上限は自分で設定できます。
	 * trueを返せば合法です。
	 * 
	 * @param kineticEnergy
	 *        Kinetic Energy.The unit is J.運動エネルギー。単位はJです。
	 * @param upperLimit
	 *        The upper limit of kinetic energy.運動エネルギーの上限です。参考値：0.98
	 * @return 
	 */
	public boolean isLegal(double kineticEnergy, double upperLimit) {
		if(kineticEnergy >= upperLimit) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * Check whether kinetic energy is legal by bullet weight and muzzle velocity. 
	 * You can set the upper limit of kinetic energy.
	 * If "true" is returned, it is legal.
	 * BB弾の重さと初速を用いて運動エネルギーが合法かどうかをチェックします。
	 * 運動エネルギーの上限は自分で設定できます。
	 * trueを返せば合法です。
	 * 
	 * @param bulletWeight
	 *        Bullet weight. The unit is g.弾丸の重さ。単位はgです。
	 * @param velocity
	 *        The speed of the bullet at the muzzle. The unit is m/s.銃口での弾丸の速度。単位はm/sです。
	 * @param upperLimit
	 *        The upper limit of kinetic energy.運動エネルギーの上限です。参考値：0.98
	 * @return
	 */
	public boolean isLegal(double bulletWeight , double velocity , double upperLimit) {
		double kineticEnergy = this.calculateKineticEnergy(bulletWeight,velocity);
		if(kineticEnergy >= upperLimit) {
			return false;
		} else {
			return true;
		}
	}
	
	/**
	 * Check whether the kinetic energy is legal. 
	 * It is suitable for Japanese law (0.98J).
	 * If "true" is returned, it is legal.
	 * 運動エネルギーが合法かどうかをチェックします。
	 * 日本の法律（0.98J）に適しております。
	 * trueを返せば合法です。
	 * 
	 * @param kineticEnergy
	 *        Kinetic Energy.The unit is J.運動エネルギー。単位はJです。
	 * @return
	 */
	public boolean isLegalJapan(double kineticEnergy) {
		if(kineticEnergy >= 0.98) {
			return false;
		} else {
			return true;
		}
	}	
	
	/**
	 * Check whether kinetic energy is legal by bullet weight and muzzle velocity. 
	 * It is suitable for Japanese law (0.98J).
	 * If "true" is returned, it is legal.
	 * BB弾の重さと初速を用いて運動エネルギーが合法かどうかをチェックします。
	 * 日本の法律（0.98J）に適しております。
	 * trueを返せば合法です。
	 * 
	 * @param bulletWeight
	 *        Bullet weight. The unit is g.弾丸の重さ。単位はgです。
	 * @param velocity
	 *        The speed of the bullet at the muzzle. The unit is m/s.銃口での弾丸の速度。単位はm/sです。
	 * @return
	 */
	public boolean isLegalJapan(double bulletWeight , double velocity) {
		double kineticEnergy = this.calculateKineticEnergy(bulletWeight,velocity);
		if(kineticEnergy >= 0.98) {
			return false;
		} else {
			return true;
		}
	}
	
	/**
	 * Check whether the specific kinetic energy (Chinese standard) is legal. 
	 * It is suitable for Chinese law (1.8J/cm^2).
	 * If "true" is returned, it is legal.
	 * 检查比动能是否合法。
	 * 适用于中国法律规定的标准（1.8J/cm^2）。
	 * 返回true则为合法。
	 * 
	 * @param specificKineticEnergy
	 *        比动能
	 * @return
	 */
	public boolean isLegalChina(double specificKineticEnergy) {
		if(specificKineticEnergy >= 1.8) {
			return false;
		} else {
			return true;
		}
	}
	
	/**
	 * Check whether the specific kinetic energy (Chinese standard) is legal by bullet weight, muzzle velocity and bullet diameter.
	 * It is suitable for Chinese law (1.8J/cm^2).
	 * Reference value of bullet diameter: 7.5(mm)
	 * If "true" is returned, it is legal.
	 * 使用子弹重量、初速、子弹直径（推荐值为7.5mm）来检查比动能是否合法。
	 * 适用于中国法律规定的标准（1.8J/cm^2）。
	 * 返回true则为合法。
	 * 
	 * @param bulletWeight
	 *        子弹重量。单位是g。
	 * @param velocity
	 *        枪口初速。单位是m/s。
	 * @param bulletDiameter
	 *        子弹直径。单位是mm。
	 * @return
	 */
	public boolean isLegalChina(double bulletWeight , double velocity,  double bulletDiameter) {
		double specificKineticEnergy = this.calculateSpecificKineticEnergyChina(bulletWeight,velocity,bulletDiameter);
		if(specificKineticEnergy >= 1.8) {
			return false;
		} else {
			return true;
		}
	}
	
	/**
	 * Calculate how long the bullet stays in the inner barrel.
	 * 弾丸がインナーバレルに留まる時間を計算します。
	 * 
	 * @param innerBarrelLength 
	 *        Inner barrel length. The unit is mm.インナーバレルの長さ。単位はmmです。
	 * @param bulletWeight 
	 *        Bullet weight. The unit is g.弾丸の重さ。単位はgです。
	 *        BB弹重量，单位克
	 * @param kineticEnergy
	 *        Kinetic Energy.The unit is J.運動エネルギー。単位はJです。
	 * @return
	 *        Time for the bullet to stay in the inner barrel. The unit is s.弾丸がインナーバレルに留まる時間、単位はsです。
	 */
	public double timeInInnerbarrel(double innerBarrelLength , double bulletWeight,  double kineticEnergy) {
		return 2*(innerBarrelLength*0.001)/Math.sqrt(2*kineticEnergy/(bulletWeight*0.001));
	}
	
	
	/**
	 * 弾道計算
	 * シキノさんの研究結果を使用（https://slpr.sakura.ne.jp/qp）
	 * 
	 * @param bulletWeight
	 * @param kineticEnergy
	 * @return
	 */
	public BallisticData ballisticCalculation(double bulletWeight, double kineticEnergy) {
		
		double[][] ballistic = new double[5000][14];
		double[][] usefulInformation = new double[2][4]; //出力結果
		double[][] landingPosition = new double[30][15]; //着弾位置記録用
		double[] initialState = new double[3]; //初期状態
		int landingNum = 0; //着弾回数
		double bulletDiameter = 5.95;
		double hopRotationAmount = 210.1; //ホップ回転量
		//double g = 9.80665;
		double a = 0.0;
		
		double humidity = 60;
		double pi = 3.14159265358979;
		int n = 12; //Number of 1st-order ODEs
		double radius = 0.001 * bulletDiameter / 2;
		double mass = 0.001 * bulletWeight;
		double energy = kineticEnergy;
		double vy = 0; //横方向の速度
		double theta = -1 * 0; //射出角度
		double height = 1.5;
		
		//運動状態
		double[] x = new double[n];
		double[] xold = new double[n];
		x[0] = 0;
		x[1] = 0;
		x[2] = height;
		x[3] = Math.sqrt(((2 * energy) / mass) - vy * vy) * Math.cos(pi * theta / 180);
		x[4] = vy;
		x[5] = Math.sqrt(((2 * energy) / mass) - vy * vy) * Math.sin(pi * theta / 180);
		x[6] = (-2) * pi * 0;
		initialState[0] = x[6];
		x[7] = -2 * pi * hopRotationAmount;
		initialState[1] = x[7];
		x[8] = -2 * pi * 0;
		initialState[2] = x[8];
		x[9] = 0;
		x[10] = 0;
		x[11] = 0;
		
		//double gravity = g;
		double temperature = 20;
		double moisture = 0.01 * humidity;
		double pressure = 101325;
		double landdist = 40; //着弾距離	
		//String wdecay = "yes"; //回転の減衰
		double dt = 0.001; //RK刻み幅
		long ntstep = 20;
		long ntmax = ballistic.length;
		
		double airDensity = rho(temperature, moisture, pressure);
		double airViscosity = eta(temperature);
		
		long irow = 0;
		double key = 0.0;
		double t = 0.0;
		ballistic[(int)irow][0] = t;
		ballistic[(int)irow][1] = x[0];
		ballistic[(int)irow][2] = x[1];
		ballistic[(int)irow][3] = x[2];
		ballistic[(int)irow][4] = x[3];
		ballistic[(int)irow][5] = x[4];
		ballistic[(int)irow][6] = x[5];
		ballistic[(int)irow][7] = x[6] / (2 * pi);
		ballistic[(int)irow][8] = (-1) * x[7] / (2 * pi);
		ballistic[(int)irow][9] = x[8] / (2 * pi);
		ballistic[(int)irow][10] = x[9];
		ballistic[(int)irow][11] = x[10];
		ballistic[(int)irow][12] = x[11];
		double ene = 0.5 * mass * (Math.pow(x[3],2) + Math.pow(x[4],2) + Math.pow(x[5],2));
		double eneold = 0.0;
		ballistic[(int)irow][13] = ene;
		double told = t;
		irow = irow + 1;
		
		for(int i=0; i<ntmax;i++) {
			rK4(n, t, dt, x, mass, radius, airDensity, airViscosity);
			if(i >= 1 && ((x[2] - height) * (xold[2] - height)) < 0) {
				if(Math.abs(x[2] - height) < Math.abs(xold[2] - height)) {
					usefulInformation[0][0] = t;
					usefulInformation[0][1] = x[0];
					usefulInformation[0][2] = Math.sqrt(Math.pow(x[3],2) + Math.pow(x[4],2) + Math.pow(x[5],2));
					ene = 0.5 * mass * (Math.pow(x[3],2) + Math.pow(x[4],2) + Math.pow(x[5],2));
					usefulInformation[0][3] = ene;
				} else {
					usefulInformation[0][0] = told;
					usefulInformation[0][1] = xold[0];
					usefulInformation[0][2] = Math.sqrt(Math.pow(xold[3],2) + Math.pow(xold[4],2) + Math.pow(xold[5],2));
					ene = 0.5 * mass * (Math.pow(xold[3],2) + Math.pow(xold[4],2) + Math.pow(xold[5],2));
					usefulInformation[0][3] = ene;
				}
			}
			
			if(i >= 1 && ((x[0] - landdist) * (xold[0] - landdist)) < 0) {
				if(key == 0) {
					//rowmax = Worksheets(sName).Cells(Rows.Count, 28).End(xlUp).Row
					landingPosition[landingNum][0] = landingNum + 1; //何回目の着弾
					
					a = (landdist - xold[0]) / (x[0] - xold[0]);
					landingPosition[landingNum][1] = told + (t - told) * a;
					landingPosition[landingNum][2] = xold[0] + (x[0] - xold[0]) * a;
					landingPosition[landingNum][3] = xold[1] + (x[1] - xold[1]) * a;
					landingPosition[landingNum][4] = xold[2] + (x[2] - xold[2]) * a;
					ene = 0.5 * mass * (Math.pow(x[3],2) + Math.pow(x[4],2) + Math.pow(x[5],2));
					eneold = 0.5 * mass * (Math.pow(xold[3],2) + Math.pow(xold[4],2) + Math.pow(xold[5],2));
					landingPosition[landingNum][5] = eneold + (ene - eneold) * a;
					
					//Initial conditions memo
					landingPosition[landingNum][6] = 2 * radius * 1000;
					landingPosition[landingNum][7] = mass * 1000;
					landingPosition[landingNum][8] = energy;
					landingPosition[landingNum][9] = height;
					landingPosition[landingNum][10] = vy;
					landingPosition[landingNum][11] = theta;
					landingPosition[landingNum][12] = initialState[0];
					landingPosition[landingNum][13] = initialState[1];
					landingPosition[landingNum][14] = initialState[2];
					
					key = 1;
				}
			}
			
			if(x[2] < 0) {
				usefulInformation[1][0] = told;
				usefulInformation[1][1] = x[0];
				usefulInformation[1][2] = Math.sqrt(Math.pow(xold[3],2) + Math.pow(xold[4],2) + Math.pow(xold[5],2));
				ene = 0.5 * mass * (Math.pow(xold[3],2) + Math.pow(xold[4],2) + Math.pow(xold[5],2));
				usefulInformation[1][3] = ene;
			}
			
			if((i+1) % ntstep == 0) {
				ballistic[(int)irow][0] = t;
				ballistic[(int)irow][1] = x[0];
				ballistic[(int)irow][2] = x[1];
				ballistic[(int)irow][3] = x[2];
				ballistic[(int)irow][4] = x[3];
				ballistic[(int)irow][5] = x[4];
				ballistic[(int)irow][6] = x[5];
				ballistic[(int)irow][7] = x[6] / (2 * pi);
				ballistic[(int)irow][8] = (-1) * x[7] / (2 * pi);
				ballistic[(int)irow][9] = x[8] / (2 * pi);
				ballistic[(int)irow][10] = x[9];
				ballistic[(int)irow][11] = x[10];
				ballistic[(int)irow][12] = x[11];
				ene = 0.5 * mass * (Math.pow(x[3],2) + Math.pow(x[4],2) + Math.pow(x[5],2));
				ballistic[(int)irow][13] = ene;
				
				irow++;
				if(x[2] < 0) {
					break;
				}
			}
			
			told = t;
			for(int j=0; j<n; j++) {
				xold[j] = x[j];
			}
				
		}
		
		BallisticData bd = new BallisticData(ballistic,usefulInformation,landingPosition,initialState);
		return bd;
	}

	private void rK4(int n, double t, double h, double[] x, double mass, double radius, double airDensity, double airViscosity) {
		double[][] k = new double[n][4];
		double[] f = new double[n];
		double[] tmp = new double[n];
		double[] c = new double[3];
		
		c[0] = 0.5;
		c[1] = 0.5;
		c[2] = 1;
		
		grk(n, t, x, f, mass, radius, airDensity, airViscosity);
		for(int i=0; i<n; i++) {
			k[i][0] = h * f[i];
		}
		for(int j=0; j<=2; j++) {
			for(int i=0; i<n; i++) {
				tmp[i] = x[i] + k[i][j] * c[j];
			}
			double tx = t + c[j] * h;
			grk(n, tx, tmp, f, mass, radius, airDensity, airViscosity);
			for(int i=0; i<n; i++) {
				k[i][j+1] = h * f[i];
			}
		}
		
		t = t + h;
		for(int i=0; i<n; i++) {
			x[i] = x[i] + (k[i][0] + k[i][3]) / 6 + (k[i][1] + k[i][2]) / 3;
		}
	}

	private void grk(int n, double t, double[] x, double[] f, double mass, double radius, double airDensity, double airViscosity) {
		double vu = 0, vd = 0, cNI = 0;
		double cm = 0;
		double pi = 3.14159265358979;
		double iw = 0.4 * mass * Math.pow(radius,2);
		double[] relv = new double[3];
		relv[0] = x[3] - x[9];
		relv[1] = x[4] - x[10];
		relv[2] = x[5] - x[11];
		
		double nv = Math.sqrt(Math.pow(relv[0],2) + Math.pow(relv[1],2) + Math.pow(relv[2],2));
		double nw = Math.sqrt(Math.pow(x[6],2) + Math.pow(x[7],2) + Math.pow(x[8],2));
		double[] l = new double[3];
		l[0] = relv[1] * x[8] - relv[2] * x[7];
		l[1] = relv[2] * x[6] - relv[0] * x[8];
		l[2] = relv[0] * x[7] - relv[1] * x[6];
		double nL = Math.sqrt(Math.pow(l[0],2) + Math.pow(l[1],2) + Math.pow(l[2],2));
		
		double cv = 0.5 * cd(nv, radius, airDensity, airViscosity) * airDensity * pi * Math.pow(radius,2) * nv;
		
		double cl = 0.12;
		
		if(nL > 0.00000000000001) {
			cm = 8 * cl * pi * Math.pow(radius,3) * airDensity * nw * nv / (3 * nL);
		} else {
			cm = 0;
		}
		
		if(nw > 0.00000000000001) {
			vu = nv * Math.sin(pi / 5.32065) - radius * nw * Math.sin(pi / 3.60475);
			vd = -nv * Math.sin(pi / 5.32065) - radius * nw * Math.sin(pi / 3.60475);
			cNI = airDensity * cf(nv, radius, airDensity, airViscosity) * pi * Math.pow(radius,3) * (Math.abs(vu) * vu + Math.abs(vd) * vd) / nw;
		} else {
			cNI = 0;
		}
		
		f[0] = x[3];
		f[1] = x[4];
		f[2] = x[5];
		
		f[3] = (-1 * cv * relv[0] - cm * l[0]) / mass;
		f[4] = (-1 * cv * relv[1] - cm * l[1]) / mass;
		f[5] = (-1 * cv * relv[2] - cm * l[2]) / mass - 9.80665;
		
		f[6] = cNI * x[6] / iw;
		f[7] = cNI * x[7] / iw;
		f[8] = cNI * x[8] / iw;
		
		f[9] = 0;
		f[10] = 0;
		f[11] = 0;
		
	}

	private double cf(double nv, double radius, double airDensity, double airViscosity) {
		double re = nv * 2 * radius * airDensity / airViscosity;
		double cf = 1.328 / Math.pow(re,0.5);
		return cf;
	}

	private double cd(double nv, double radius, double airDensity, double airViscosity) {
		double re = nv * 2 * radius * airDensity / airViscosity;
		double c1 = 24 / re;
		double  c2 = re / 5;
		c2 = 2.6 * c2 / (1 + Math.pow(c2,1.52));
		double  c3 = re / 263000;
		c3 = 0.411 * Math.pow(c3,-7.94) / (1 + Math.pow(c3,-8));
		double  c4 = Math.pow(re,0.8) / 461000;
				 
		double  cd = c1 + c2 + c3 + c4;
		return cd;
	}

	private double eta(double temperature) {
		double eta = 0.000001487 * Math.pow((temperature + 273.15),1.5) / ((temperature + 273.15) + 117);
		return eta;
	}

	private double rho(double temperature, double moisture, double pressure) {
		double et = 6.1078 * Math.pow(10,(7.5 * temperature / (temperature + 237.3)));
		et = 100 * et;
		et = moisture * et;
		double rhoair = 0.0034856447 * pressure / (temperature + 273.15 - 0.67);
		double rho = rhoair * (1 - 0.378 * et / pressure);
		return rho;
	}
}
