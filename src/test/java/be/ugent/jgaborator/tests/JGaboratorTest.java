package be.ugent.jgaborator.tests;

import be.tarsos.dsp.AudioDispatcher;
import be.tarsos.dsp.io.TarsosDSPAudioFloatConverter;
import be.tarsos.dsp.io.jvm.AudioDispatcherFactory;
import be.tarsos.dsp.io.jvm.JVMAudioInputStream;
import be.ugent.jgaborator.JGaborator;
import org.junit.jupiter.api.Test;

import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.UnsupportedAudioFileException;
import java.io.IOException;
import java.net.URL;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class JGaboratorTest {
    @Test
    public void runOnFile(){
        final float[]  audioSamples = audioBufferFile("44.1kHz_440Hz_1s.wav",44100);

        try {
            AudioDispatcher adp = AudioDispatcherFactory.fromFloatArray(audioSamples,44100,1024,0);

            final JGaborator zsazsa = new JGaborator(1024 , 44100, 24, 220,880,440.5  ,512);

            adp.addAudioProcessor(zsazsa);
            adp.run();
            List<float[]> coefficients = zsazsa.getCoefficents();
            assertTrue(coefficients.size() > 0, "Coeffients not calculated");

            //Find the max bin index in the coefficients in the middle of the audio file
            float[] center = coefficients.get(coefficients.size()/2);
            int maxIndex = 0;
            double maxValue = -11000;
            for(int i = 0 ; i < center.length ; i++){
                if(center[i] > maxValue){
                    maxIndex = i;
                    maxValue = center[i];
                }
            }
            //should be close to band center close to 440Hz
            assertEquals(440,zsazsa.bandCenters(maxIndex),0.501);
        } catch (UnsupportedAudioFileException e) {
            throw new RuntimeException(e);
        }

    }

    private static float[] audioBufferFile(String file,int lengthInSamples){
        float[] buffer = new float[lengthInSamples];
        try {
            final URL url = ClassLoader.getSystemResource(file);
            AudioInputStream audioStream = AudioSystem.getAudioInputStream(url);
            AudioFormat format = audioStream.getFormat();
            TarsosDSPAudioFloatConverter converter = TarsosDSPAudioFloatConverter.getConverter(JVMAudioInputStream.toTarsosDSPFormat(format));
            byte[] bytes = new byte[lengthInSamples * format.getSampleSizeInBits()];
            audioStream.read(bytes);
            converter.toFloatArray(bytes, buffer);
        } catch (IOException e) {
            throw new Error("Test audio file should be present.");
        } catch (UnsupportedAudioFileException e) {
            throw new Error("Test audio file format should be supported.");
        }
        return buffer;
    }
}
