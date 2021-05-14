/** 
* A calculator for Airsoft.
* This is the main window used to start the program.
* 
* @author  Karasukaigan 
* @date    2021/5/15 
* @version 1.0.0
*/

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.EventQueue;
import java.awt.Graphics;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.JDesktopPane;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JButton;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class MainWindow extends JFrame {

	private static final long serialVersionUID = 1L;
	private JPanel contentPane;
	private JTextField textField;
	private JTextField textField_1;
	public double bulletWeight = 0.0, velocity = 0.0;
	public double kineticEnergy = 0.0;
	public String isLegal = "";
	public BallisticData bd;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					MainWindow frame = new MainWindow();
					frame.setLocationRelativeTo(null);
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the frame.
	 * @throws UnsupportedLookAndFeelException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 * @throws ClassNotFoundException 
	 */
	public MainWindow() throws ClassNotFoundException, InstantiationException, IllegalAccessException, UnsupportedLookAndFeelException {
		UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		setTitle("エアガン専用電卓");
		setResizable(false);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 654, 322);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.setLayout(new BorderLayout(0, 0));
		setContentPane(contentPane);
		
		JDesktopPane desktopPane = new JDesktopPane();
		desktopPane.setBackground(UIManager.getColor("Button.background"));
		contentPane.add(desktopPane, BorderLayout.CENTER);
		
		JLabel lblNewLabel = new JLabel("運動エネルギー計算");
		lblNewLabel.setBounds(20, 20, 120, 15);
		desktopPane.add(lblNewLabel);
		
		JLabel lblNewLabel_1 = new JLabel("BB弾重さ：");
		lblNewLabel_1.setBounds(20, 46, 60, 15);
		desktopPane.add(lblNewLabel_1);
		
		textField = new JTextField();
		textField.setText("0.2");
		textField.setBounds(82, 43, 49, 21);
		desktopPane.add(textField);
		textField.setColumns(10);
		
		JLabel lblNewLabel_2 = new JLabel("g");
		lblNewLabel_2.setBounds(136, 45, 21, 15);
		desktopPane.add(lblNewLabel_2);
		
		JLabel lblNewLabel_3 = new JLabel("初速：");
		lblNewLabel_3.setBounds(156, 46, 36, 15);
		desktopPane.add(lblNewLabel_3);
		
		textField_1 = new JTextField();
		textField_1.setColumns(10);
		textField_1.setBounds(193, 43, 49, 21);
		desktopPane.add(textField_1);
		
		JLabel lblNewLabel_2_1 = new JLabel("m/s");
		lblNewLabel_2_1.setBounds(248, 46, 21, 15);
		desktopPane.add(lblNewLabel_2_1);
		
		JLabel lblNewLabel_4 = new JLabel("運動エネルギー：");
		lblNewLabel_4.setBounds(388, 46, 240, 15);
		desktopPane.add(lblNewLabel_4);
		
		JLabel lblNewLabel_5 = new JLabel("弾道計算");
		lblNewLabel_5.setBounds(20, 79, 54, 15);
		desktopPane.add(lblNewLabel_5);
		
		JLabel lblNewLabel_6 = new JLabel("射出高さに着弾する距離：0m");
		lblNewLabel_6.setBounds(20, 103, 212, 15);
		desktopPane.add(lblNewLabel_6);
		
		
		
		JButton btnNewButton = new JButton("計算");
		btnNewButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		btnNewButton.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				AirsoftCalculator ac = new AirsoftCalculator();
				bulletWeight = Double.parseDouble(textField.getText());
				velocity = Double.parseDouble(textField_1.getText());
				kineticEnergy = ac.calculateKineticEnergyRounding(bulletWeight,velocity);
				if(ac.isLegal(kineticEnergy, 0.989) == true) {
					isLegal = "合法";
				} else {
					isLegal = "違法";
				}
				lblNewLabel_4.setText("運動エネルギー：" + kineticEnergy + "J （" + isLegal + "）");
				
				bd = ac.ballisticCalculation(bulletWeight, kineticEnergy);
				lblNewLabel_6.setText("射出高さに着弾する距離："+ ((double)Math.round(bd.usefulInformation[0][1] * 1000) / 1000) +"m");
				
				JPanel panel = new JPanel() {
					private static final long serialVersionUID = 1L;
					
					public void paint(Graphics graphics) {
						super.paint(graphics);
						graphics.setColor(Color.blue);
						
						double x1 = 0, y1 = 0;
						double x2 = 0, y2 = 0;
		                
						for(int i=0; i<bd.ballistic.length-1; i++) {
							if(bd.ballistic[i][1] % 10 <= 0.6) {
								graphics.drawString((int)bd.ballistic[i][1]+"m", (int)(bd.ballistic[i][1] / 45 * this.getWidth()), (int)(bd.ballistic[i][3] / 2 * this.getHeight() * (-1) + 200 - 10));
							}
							if(bd.ballistic[i][3] >= bd.ballistic[0][3] - 0.1) {
								graphics.setColor(Color.blue);
								x1 = bd.ballistic[i][1] / 45 * this.getWidth();
								y1 = bd.ballistic[i][3] / 2 * this.getHeight() * (-1) + 200;
								x2 = bd.ballistic[i+1][1] / 45 * this.getWidth();
								y2 = bd.ballistic[i+1][3] / 2 * this.getHeight() * (-1) + 200;
								graphics.drawLine((int)x1, (int)y1, (int)x2, (int)y2);
							} else {
								graphics.setColor(Color.red);
								x1 = bd.ballistic[i][1] / 45 * this.getWidth();
								y1 = bd.ballistic[i][3] / 2 * this.getHeight() * (-1) + 200;
								x2 = bd.ballistic[i+1][1] / 45 * this.getWidth();
								y2 = bd.ballistic[i+1][3] / 2 * this.getHeight() * (-1) + 200;
								graphics.drawLine((int)x1, (int)y1, (int)x2, (int)y2);
							}
							
							System.out.println((int)(bd.ballistic[i][1]/50*579)+ ","+(int)(bd.ballistic[i][3]/1.5*192));
							
						}
						
					}
				};
				panel.setBounds(20, 128, 597, 192);
				desktopPane.add(panel);
			}
		});
		btnNewButton.setBounds(283, 42, 95, 23);
		desktopPane.add(btnNewButton);	
		
	}
}
