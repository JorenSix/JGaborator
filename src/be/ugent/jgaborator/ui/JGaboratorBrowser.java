package be.ugent.jgaborator.ui;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.Random;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import be.tarsos.dsp.AudioDispatcher;
import be.tarsos.dsp.io.jvm.AudioDispatcherFactory;
import be.tarsos.dsp.ui.AxisUnit;
import be.tarsos.dsp.ui.CoordinateSystem;
import be.tarsos.dsp.ui.LinkedPanel;
import be.tarsos.dsp.ui.ViewPort;
import be.tarsos.dsp.ui.ViewPort.ViewPortChangedListener;
import be.tarsos.dsp.ui.layers.BackgroundLayer;
import be.tarsos.dsp.ui.layers.DragMouseListenerLayer;
import be.tarsos.dsp.ui.layers.SelectionLayer;
import be.tarsos.dsp.ui.layers.TimeAxisLayer;
import be.tarsos.dsp.ui.layers.ZoomMouseListenerLayer;
import be.ugent.jgaborator.JGaborator;


public class JGaboratorBrowser  extends JFrame{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1973676584543780097L;

	private final JPanel fingerprintPanel;

	private String path;
	
	private double minFrequency,maxFrequency,refFrequency; //in Hz
	private int bandsPerOctave;
	private int sampleRate; // in Hz
	private int resolution; //in audio samples
	private int stepSize;//audio block step size in samples
	
	
	public JGaboratorBrowser() {
		
		minFrequency = 110;
		maxFrequency = 3520;
		refFrequency = 440;
		bandsPerOctave = 12;
		sampleRate = 8000;
		resolution = 64; 
		stepSize = 4096* 2 * 2;
		
		this.setLayout(new BorderLayout());
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		this.setTitle("Gaborator Visualizer");
		
		fingerprintPanel = new JPanel();
		fingerprintPanel.setLayout(new GridLayout(0,1));
		fingerprintPanel.add(emptyFeaturePanel());
		
		new FileDrop(null, fingerprintPanel, /*dragBorder,*/ new FileDrop.Listener(){   
			public void filesDropped( java.io.File[] files ){   
				for( int i = 0; i < files.length; i++) {   
					final File fileToAdd = files[i];
					new Thread(new Runnable(){
						@Override
						public void run() {					
		                	addAudio(fileToAdd.getAbsolutePath());
						}}).start();
                }
			}
        });
		
		this.add(fingerprintPanel,BorderLayout.CENTER);
		this.add(controlComponent(),BorderLayout.WEST);
	}
	

	
	public Component controlComponent() {
		
		JPanel container = new JPanel(new BorderLayout(10,10));
		container.setBorder(new EmptyBorder(10, 10, 10, 10));
		JPanel p = new JPanel(new GridLayout(0, 2, 10, 10));
		
		JLabel blockSizeLabel = new JLabel("Audio block size");
		blockSizeLabel.setToolTipText("in audio samples");
		JSpinner blockSizeSpinner = new JSpinner(new SpinnerNumberModel(stepSize,1,2000000,1));
		blockSizeSpinner.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				stepSize =  ((Integer) ((JSpinner) e.getSource()).getValue()).intValue();
			}
		});
		
		p.add(blockSizeLabel);
		p.add(blockSizeSpinner);
		
		
		JLabel resolutionLabel = new JLabel("Time bin resolution");
		resolutionLabel.setToolTipText("in audio samples");
		JSpinner resolutionSpinner = new JSpinner(new SpinnerNumberModel(resolution,1,4096*2*2,1));
		resolutionSpinner.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				resolution =  ((Integer) ((JSpinner) e.getSource()).getValue()).intValue();
			}
		});
		p.add(resolutionLabel);
		p.add(resolutionSpinner);
		
		
		JLabel samplerateLabel = new JLabel("Audio sample rate");
		samplerateLabel.setToolTipText("in Hertz");
		JSpinner samplerateSpinner = new JSpinner(new SpinnerNumberModel(sampleRate,1,96000,1));
		samplerateSpinner.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				sampleRate =  ((Integer) ((JSpinner) e.getSource()).getValue()).intValue();
			}
		});
		p.add(samplerateLabel);
		p.add(samplerateSpinner);
		
		
		JLabel bandsPerOctaveLabel = new JLabel("Bands per octave");
		JSpinner bandsPerOctaveSpinner = new JSpinner(new SpinnerNumberModel(bandsPerOctave,1,120,1));
		bandsPerOctaveSpinner.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				bandsPerOctave =  ((Integer) ((JSpinner) e.getSource()).getValue()).intValue();
			}
		});
		p.add(bandsPerOctaveLabel);
		p.add(bandsPerOctaveSpinner);
		
		JLabel minFrequencyLabel = new JLabel("Min frequency");
		minFrequencyLabel.setToolTipText("in Hertz");
		JSpinner minFrequencySpinner = new JSpinner(new SpinnerNumberModel(minFrequency,1,25000,1));
		minFrequencySpinner.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				minFrequency =  ((Double) ((JSpinner) e.getSource()).getValue()).intValue();
			}
		});		
		p.add(minFrequencyLabel);
		p.add(minFrequencySpinner);
		
		JLabel maxFrequencyLabel = new JLabel("Max frequency");
		maxFrequencyLabel.setToolTipText("in Hertz");
		JSpinner maxFrequencySpinner = new JSpinner(new SpinnerNumberModel(maxFrequency,1,25000,1));
		maxFrequencySpinner.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				maxFrequency = ((Double) ((JSpinner) e.getSource()).getValue()).intValue();
			}
		});		
		p.add(maxFrequencyLabel);
		p.add(maxFrequencySpinner);
		
		JLabel refFrequencyLabel = new JLabel("Ref frequency");
		refFrequencyLabel.setToolTipText("in Hertz");
		
		JSpinner refFrequencySpinner = new JSpinner(new SpinnerNumberModel(refFrequency,1,25000,1));
		refFrequencySpinner.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				refFrequency = ((Double) ((JSpinner) e.getSource()).getValue()).intValue();
			}
		});		
		p.add(refFrequencyLabel);
		p.add(refFrequencySpinner);	
		
		
		JButton clearButton = new JButton("Clear");
		clearButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				Runnable uiRunnable = new Runnable() {
					@Override
					public void run() {
						fingerprintPanel.removeAll();
						fingerprintPanel.add(emptyFeaturePanel());
						fingerprintPanel.validate();
					}
				};				
				// run ui stuff on ui thread.
				SwingUtilities.invokeLater(uiRunnable);
			}
		});
		
		
		JButton recalcButton = new JButton("Recalc");
		recalcButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				addAudio(path);
			}
		});
		
		p.add(clearButton);
		p.add(recalcButton);
		
		container.add(p,BorderLayout.NORTH);
		container.add(new JPanel(),BorderLayout.CENTER);
		
		return container;
	}
	
	private void addAudio(final String path){
		this.path = path;
		
		Runnable clearRunnable = new Runnable() {
			@Override
			public void run() {
				fingerprintPanel.removeAll();
				fingerprintPanel.add(emptyFeaturePanel());
				fingerprintPanel.validate();
				//Validating a container means laying out its subcomponents:
			}
		};
		// run ui stuff on ui thread.
		SwingUtilities.invokeLater(clearRunnable);
		
		//audio processing on a new thread
		new Thread(new Runnable() {

			@Override
			public void run() {
				final JGaborator zsazsa = new JGaborator(stepSize, sampleRate, bandsPerOctave, minFrequency,maxFrequency,refFrequency,resolution);
				AudioDispatcher ad = AudioDispatcherFactory.fromPipe(path, sampleRate, stepSize, 0);
				ad.addAudioProcessor(zsazsa);
				ad.run();
				
				System.out.println("Finished analysis of " + path);
				
				//after audio processing, update ui (on the ui thread)
				Runnable clearRunnable = new Runnable() {
					@Override
					public void run() {
						fingerprintPanel.removeAll();
						Component featurePanel = createFeaturePanel(zsazsa);
						fingerprintPanel.add(featurePanel);
						fingerprintPanel.validate();
						//Validating a container means laying out its subcomponents:
					}
				};
				// run ui stuff on ui thread.
				SwingUtilities.invokeLater(clearRunnable);
			}},"Audio processing").start();;		
	}
	
	private Component emptyFeaturePanel(){
		final CoordinateSystem cs = new CoordinateSystem(AxisUnit.FREQUENCY, 3500, 11900);
		final LinkedPanel frequencyDomainPanel = new LinkedPanel(cs);
		frequencyDomainPanel.getViewPort().addViewPortChangedListener(new ViewPortChangedListener() {
			@Override
			public void viewPortChanged(ViewPort newViewPort) {
				frequencyDomainPanel.repaint();
			}
		});	
		frequencyDomainPanel.addLayer(new ZoomMouseListenerLayer());
		frequencyDomainPanel.addLayer(new DragMouseListenerLayer(cs));
		frequencyDomainPanel.addLayer(new BackgroundLayer(cs));
		frequencyDomainPanel.addLayer(new FrequencyAxisLayer(cs));
		frequencyDomainPanel.addLayer(new TimeAxisLayer(cs));
		frequencyDomainPanel.addLayer(new SelectionLayer(cs));
		return frequencyDomainPanel;
	}
	
	private Component createFeaturePanel(JGaborator zsazsa) {
		
		final CoordinateSystem cs = new CoordinateSystem(AxisUnit.FREQUENCY, 3500, 11900);
		final LinkedPanel frequencyDomainPanel = new LinkedPanel(cs);
		frequencyDomainPanel.getViewPort().addViewPortChangedListener(new ViewPortChangedListener() {
			@Override
			public void viewPortChanged(ViewPort newViewPort) {
				frequencyDomainPanel.repaint();
			}
		});
		
		frequencyDomainPanel.addLayer(new ZoomMouseListenerLayer());
		frequencyDomainPanel.addLayer(new DragMouseListenerLayer(cs));
		frequencyDomainPanel.addLayer(new BackgroundLayer(cs));
		frequencyDomainPanel.addLayer(new GaborLayer(cs,zsazsa));
		frequencyDomainPanel.addLayer(new FrequencyAxisLayer(cs));
		frequencyDomainPanel.addLayer(new TimeAxisLayer(cs));
		frequencyDomainPanel.addLayer(new SelectionLayer(cs));
		
		return frequencyDomainPanel;
	}
	
	
	public static void main(final String[] args) {
		if (args.length == 0) {
			JFrame frame = null;
			frame = new JGaboratorBrowser();
			frame.pack();
			frame.setSize(800,550);
			frame.setLocationRelativeTo(null);
			frame.setVisible(true);
		}else {
			for(int i = 0 ; i < 2 ; i++) {
				
				int index = 0;
				final boolean singleThreaded = false;
				for(final String path : args ) {
					if(!new File(path).exists())
						continue;
					
					index++;
					final int count = index;
					final int in = i+1;
					
					//start all threads simultaneously
					Thread t = new Thread(new Runnable() {
						@Override
						public void run() {
							if(!singleThreaded) {
								try {
									Thread.sleep((long) (new Random().nextFloat() * 10));
								} catch (InterruptedException e) {
									e.printStackTrace();
								}
							}
							
							double minFrequency = 110;
							double maxFrequency = 3520;
							double refFrequency = 440;
							int bandsPerOctave = 12;
							double sampleRate = 8000;
							int resolution = 64; 
							int stepSize = 4096* 2 * 2;
							
							System.out.printf("%d %d/%d START %s\n",in, count,args.length,path);
							
							final JGaborator zsazsa = new JGaborator(stepSize, sampleRate, bandsPerOctave, minFrequency,maxFrequency,refFrequency,resolution);
							AudioDispatcher ad = AudioDispatcherFactory.fromPipe(path, (int)sampleRate, stepSize, 0);
							ad.addAudioProcessor(zsazsa);
							ad.run();
							
							System.out.printf("%d %d/%d STOP %s\n",in,count,args.length,path);
						}
					});
					if(singleThreaded)
						t.run();
					else
						t.start();
				}
			}
		}
	}

}
