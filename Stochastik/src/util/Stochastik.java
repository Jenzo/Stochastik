package util;

import java.math.BigInteger;
import java.text.DecimalFormat;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;

public class Stochastik {


	public static BigInteger fakultaet(final long n) {

		return n == 0 ? BigInteger.valueOf(1)
				: fakultaet(BigInteger.valueOf(n - 1).longValueExact()).multiply(BigInteger.valueOf(n));
	}

	public static BigInteger nUeberK(final long n, final long k) {

		return fakultaet(n).divide((fakultaet(k).multiply(fakultaet(n - k))));
	}

	public static double hyperGeometrischeVeteilung(long r, long s, long n, long k) {

		if (k > n) {
			throw new IllegalArgumentException("k darf nicht größer als n sein");
		}

		if (k > r) {
			throw new IllegalArgumentException("k darf nicht größer als r oder s sein");
		}

		if (n > r + s) {
			throw new IllegalArgumentException("n darf nicht größer als r + s sein");
		}

		long gesamt = r + s;

		double wahrscheinlichkeit_Ak = nUeberK(r, k).multiply(nUeberK(s, n - k)).doubleValue();
		double wahrscheinlichkeit_alle_Ergebnisse = nUeberK(gesamt, n).doubleValue();

		return wahrscheinlichkeit_Ak / wahrscheinlichkeit_alle_Ergebnisse;

	}

	public static double binominalVerteilung(long r, long s, long n, long k) {

		if (k > n) {
			throw new IllegalArgumentException("k darf nicht größer als n sein");
		}

		long kVonN_roteKugeln = nUeberK(n, k).longValueExact();
		double anteilRot = Math.pow(((r / r + s)), k);
		double anteilSchwarz = Math.pow(((s / s + r)), n - k);

		double ak = kVonN_roteKugeln * anteilRot * anteilSchwarz;
		double pOmega = Math.pow((r + s), n);

		return ak / pOmega;
	}

	public static double bedingteWahrscheinlichkeit(final long gesamt, final long anteil, final long zuege,
			final long erwartung) {

		final double omega = nUeberK(gesamt, zuege).longValueExact();
		final double ereignis = nUeberK(gesamt - anteil, zuege).multiply(nUeberK(anteil, erwartung)).longValueExact();

		return ereignis / omega;
	}

	public static double BernulliTuple(int gesamt, double wahrscheinlichkeit, int anteil) {

		return nUeberK(gesamt, anteil).doubleValue() * Math.pow(wahrscheinlichkeit, anteil)
				* Math.pow(1 - wahrscheinlichkeit, gesamt - anteil);
	}

	public static double BernulliAbschätzung(int gesamt, double wahrscheinlichkeit, int anteil) {
		double sigma = gesamt * wahrscheinlichkeit * (1 - wahrscheinlichkeit);
		double mü = gesamt * wahrscheinlichkeit;

		return (1 / Math.sqrt(sigma)) * gammaVonX((anteil - mü) / Math.sqrt(sigma));

	}

	private static double gammaVonX(double x) {
		return (1 / (Math.sqrt(2 * Math.PI))) * Math.pow(Math.E, -(x * x) / 2);
	}

	public static double PoissonNährung(BernoulliTuple b) {
		double my = b.getN() * b.getP();
		long fakK = fakultaet(b.getK()).longValueExact();

		return (Math.pow(my, b.getK()) / fakK) * Math.pow(Math.E, -my);

	}

	public static double PoissonNährung(int n, double p, int max) {
		double result = 0;
		for (int k = 0; k < max; k++) {
			result += PoissonNährung(new BernoulliTuple(n, p, k));
		}
		return result;
	}

	public static double x(int n, double p, int max, int min) {
		double µ = n * p;
		double sigma = µ * (1 - p);

		double k2 = normalverteilung(max, µ, sigma);
		double k1 = normalverteilung(min, µ, sigma);

		return bla(k2) - bla(k1);
	}

	public static double bla(double x) {
		System.out.println(Math.pow(Math.E, Math.pow(x, 2)));
		return (1 / Math.sqrt(2 * Math.PI)) * Math.pow(Math.E, Math.pow(-0.5 * x, 2));
	}

	private static double normalverteilung(int k, double µ, double sigma) {
		return (k - µ + 0.5) / Math.sqrt(sigma);
	}

	public static class BernoulliTuple {

		private final int n;
		private final double p;
		private final int k;

		public BernoulliTuple(int n, double p, int k) {
			this.n = n;
			this.p = p;
			this.k = k;
		}

		public int getN() {
			return n;
		}

		public double getP() {
			return p;
		}

		public int getK() {
			return k;
		}
	}

	public static double GaußscheNormalverteilungTabellenwert(double x) {

		final double relativeAccuracy = 1e-4;
		final double absoluteAccuracy = 1e-8;
		final int minimalIterationCount = 3;
		final int maximalIterationCount = 32;

		UnivariateIntegrator integrator = new SimpsonIntegrator(relativeAccuracy, absoluteAccuracy,
				minimalIterationCount, maximalIterationCount);
		UnivariateFunction f = new UnivariateFunction() {

			@Override
			public double value(double x) {
				return Math.exp(-Math.pow(x, 2) / 2);
			}
		};

		return (1 / Math.sqrt(Math.PI * 2)) * integrator.integrate(10000, f, -500, x);
	}

	public static String getProportionaleFunktion(Punkt[] punkte) {
		double m;
		double sumXiYi = 0, sumXipow = 0;

		for (int i = 0; i < punkte.length; i++) {
			sumXiYi += punkte[i].x * punkte[i].y;
			sumXipow += Math.pow(punkte[i].x, 2);
		}

		m = sumXiYi / sumXipow;

		return "f(x) = " + m + "*x";

	}

	public static String getLineareFuntion(Punkt[] punkte) {

		double a, b;
		int n = punkte.length;

		double xQuer = 0, yQuer = 0;
		for (int i = 0; i < n; i++) {
			xQuer += punkte[i].x / n;
			yQuer += punkte[i].y / n;
		}

		double sumXiYi = 0, sumXiPow = 0;
		for (int i = 0; i < n; i++) {
			sumXiYi += punkte[i].x * punkte[i].y;
			sumXiPow += Math.pow(punkte[i].x, 2);
		}

		a = (sumXiYi - n * xQuer * yQuer) / (sumXiPow - n*Math.pow(xQuer, 2));
		b = yQuer - a * xQuer;

		return "f(x) = " + a + "*x" + " + " + b;
	}
	
	public static String getExponentialFunktion(Punkt[] punkte) {
		
		DecimalFormat format = new DecimalFormat("#.####");
		int n = punkte.length;
		Punkt[] neuePunkte = new Punkt[n];
		for(int i = 0 ; i < n ; i++) {
			neuePunkte[i] = new Punkt(punkte[i].x, Math.log(punkte[i].y));
		}
		
		double xQuer = 0 , yQuer = 0;
		double sumXiYi = 0 , sumXiPow = 0 ;
		for(int i = 0 ; i < n ; i++) {
			xQuer += neuePunkte[i].x / n;
			yQuer += neuePunkte[i].y / n;
			sumXiYi += neuePunkte[i].x * neuePunkte[i].y;
			sumXiPow += Math.pow(neuePunkte[i].x, 2);
			
		}
		
		Punkt mittelpunkt = new Punkt(xQuer, yQuer);
		double a, b;
		a = (sumXiYi - (n * xQuer * yQuer)) / (sumXiPow - (n * Math.pow(xQuer, 2)));
		b = yQuer - (a * xQuer);
		
		double neuesA = Math.exp(b);
		double neuesB = a;
		System.out.println(String.format("Mittelpunkt : (%s, %s)", mittelpunkt.x, format.format(mittelpunkt.y)));
		return "f(x) = a*e^bx = " + neuesA + "*e^" + neuesB + "*x";
		
	}
	
	public static String getPotenzFunktion(Punkt[] punkte) {
		DecimalFormat format = new DecimalFormat("#.####");
		int n = punkte.length;
		Punkt[] neuePunkte = new Punkt[n];
		for(int i = 0 ; i < n ; i++) {
			neuePunkte[i] = new Punkt(Math.log(punkte[i].x), Math.log(punkte[i].y));
		}
		
		double xQuer = 0 , yQuer = 0;
		double sumXiYi = 0 , sumXiPow = 0 ;
		for(int i = 0 ; i < n ; i++) {
			xQuer += neuePunkte[i].x / n;
			yQuer += neuePunkte[i].y / n;
			sumXiYi += neuePunkte[i].x * neuePunkte[i].y;
			sumXiPow += Math.pow(neuePunkte[i].x, 2);
			
		}
		
		Punkt mittelpunkt = new Punkt(xQuer, yQuer);
		double a, b;
		a = (sumXiYi - (n * xQuer * yQuer)) / (sumXiPow - (n * Math.pow(xQuer, 2)));
		b = yQuer - (a * xQuer);
		
		double neuesA = Math.exp(b);
		double neuesB = a;
		System.out.println(String.format("Mittelpunkt : (%s, %s)", mittelpunkt.x, format.format(mittelpunkt.y)));
		return "f(x) = a*x^b = " + neuesA + "*x^" + neuesB;
		
	}

}
