package net.sytes.botg.array.math;

public class Vec2Mat {

	// Suppress default constructor for noninstantiability
	private Vec2Mat() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}	
	
	public static double[][] spectrogram(double[] x, double sampleRate, int fftSampleSize, int overlapFactor) {

		int numSamples = x.length;

		double[][] windows = overlapWindows(x, fftSampleSize, overlapFactor);
		
		int w = windows.length;
		
		/*
		
		numFrames = numSamples / fftSampleSize;
		framesPerSecond = (int) (numFrames / wave.getLengthInSeconds());

		// set signals for fft
		WindowFunction window = new WindowFunction();
		double[] win = window.generate(WindowType.HAMMING, fftSampleSize);

		double[][] signals = new double[numFrames][];
		for (int f = 0; f < numFrames; f++) {
			signals[f] = new double[fftSampleSize];
			int startSample = f * fftSampleSize;
			for (int n = 0; n < fftSampleSize; n++) {
				signals[f][n] = amplitudes[startSample + n] * win[n];
			}
		}

		absoluteSpectrogram = new double[numFrames][];
		// for each frame in signals, do fft on it
		FastFourierTransform fft = new FastFourierTransform();
		for (int i = 0; i < numFrames; i++) {
			absoluteSpectrogram[i] = fft.getMagnitudes(signals[i]);
		}

		if (absoluteSpectrogram.length > 0) {

			numFrequencyUnit = absoluteSpectrogram[0].length;
			frequencyBinSize = (double) wave.getWaveHeader().getSampleRate() / 2 / numFrequencyUnit; // frequency could
																										// be caught
																										// within the
																										// half of
																										// nSamples
																										// according to
																										// Nyquist
																										// theory
			frequencyRange = wave.getWaveHeader().getSampleRate() / 2;

			// normalization of absoultSpectrogram
			spectrogram = new double[numFrames][numFrequencyUnit];

			// set max and min amplitudes
			double maxAmp = Double.MIN_VALUE;
			double minAmp = Double.MAX_VALUE;
			for (int i = 0; i < numFrames; i++) {
				for (int j = 0; j < numFrequencyUnit; j++) {
					if (absoluteSpectrogram[i][j] > maxAmp) {
						maxAmp = absoluteSpectrogram[i][j];
					} else if (absoluteSpectrogram[i][j] < minAmp) {
						minAmp = absoluteSpectrogram[i][j];
					}
				}
			}

			// normalization
			// avoiding divided by zero
			double minValidAmp = 0.00000000001F;
			if (minAmp == 0) {
				minAmp = minValidAmp;
			}

			double diff = Math.log10(maxAmp / minAmp); // perceptual difference
			for (int i = 0; i < numFrames; i++) {
				for (int j = 0; j < numFrequencyUnit; j++) {
					if (absoluteSpectrogram[i][j] < minValidAmp) {
						spectrogram[i][j] = 0;
					} else {
						spectrogram[i][j] = (Math.log10(absoluteSpectrogram[i][j] / minAmp)) / diff;
						// System.out.println(spectrogram[i][j]);
					}
				}
			}
		}
		 */
		return null;
	}
	
	/**
	 * separates the given double array {@code ar}
	 * @param ar
	 * @param windowSize
	 * @param overlap
	 * @return
	 */
	public static double[][] overlapWindows(double[] ar, int windowSize, int overlap){
		int n = ar.length;
		if (n < windowSize) {
			throw new IllegalArgumentException("ar.length is smaller than windowSize");
		}
		if (overlap < 0 || overlap > windowSize - 1) {
			throw new IllegalArgumentException("overlap must be in range [0, windowSize - 1]");
		}
		int o;
		int wins;
		int rem;
		if (overlap == 0) {
			wins = n / windowSize;
			rem = n % windowSize;
			if (rem > 0) {
				++wins;
			}
		} else {
			wins = (n - overlap) / (windowSize - overlap);
			rem =  (n - overlap) % (windowSize - overlap);
			if (rem > 0) {
				++wins;
			}
		}
		
		double[][] windows = new double[wins][windowSize];
		int w = 0;
		int j = 0;
		for (int i = 0; i < n; i++) {
			if (j >= windowSize) {
				++w;
				j = 0;
				//i = i - (windowSize - o);
				i = i - overlap;
				if (w == wins - 1) {
					i = n - windowSize;
				}
			}
			windows[w][j] = ar[i];
			++j;
		}
		return windows;
	}
	
	/**
	 * separates the given double array into windows of size windowSize
	 * <br>if dropUneven is true, then the remaining last ar.length % windowSize elements are dropped
	 * <br>dimensions --> double[win][windowSize]
	 * @param ar
	 * @param windowSize
	 * @param dropUneven
	 * @return
	 */
	public static double[][] fixedWindows(double[] ar, int windowSize, boolean dropUneven){
		if (ar == null) {
			return null;
		}
		int win = 0;
		int rem = 0;
		if (dropUneven) {
			win = ar.length / windowSize;
		} else {
			win = ar.length / windowSize;
			rem =  ar.length % windowSize;
			if (rem > 0) {
				++win;
			}
		}
		double[][] wAr = new double[win][];
		for (int i = 0; i < win; i++) {
			double[] dAr = null;
			if (i == win - 1) {
				if (dropUneven) {
					dAr = new double[windowSize];
					System.arraycopy(ar, i * windowSize, dAr, 0, windowSize);
				} else {
					if (rem == 0) {
						dAr = new double[windowSize];
						System.arraycopy(ar, i * windowSize, dAr, 0, windowSize);
					} else {
						dAr = new double[rem];
						System.arraycopy(ar, i * windowSize, dAr, 0, rem);
					}
				}
			} else {
				dAr = new double[windowSize];
				System.arraycopy(ar, i * windowSize, dAr, 0, windowSize);
			}
			wAr[i] = dAr;
		}
		return wAr;
	}
	
}
