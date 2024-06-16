package net.sytes.botg.array.math;

public enum Feature {
	
	SUM, MIN, MAX, MEAN, SPAN, MEDIAN, RMS, RMSMEAN, VARIANCE, SKEWNESS, KURTOSIS, CREST, NORM;
	
	public static String[] names() {
		Feature[] features = Feature.values();
		String[] names = new String[features.length];
		for (int f = 0; f < names.length; f++) {
			names[f] = features[f].name();
		}
		return names;
	}
	
}
