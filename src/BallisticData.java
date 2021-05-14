
public class BallisticData {
	public double[][] ballistic = new double[5000][14]; //運動状態
	public double[][] usefulInformation = new double[2][4]; //出力結果
	public double[][] landingPosition = new double[30][15]; //着弾位置記録用
	public double[] initialState = new double[3]; //初期状態
	
	public BallisticData(double[][] ballistic, double[][] usefulInformation, double[][] landingPosition, double[] initialState) {
		this.ballistic = ballistic;
		this.usefulInformation = usefulInformation;
		this.landingPosition = landingPosition;
		this.initialState = initialState;
	}
}
