RenyiEntropy : MultiOutUGen {

	*kr {
		arg fft, fftsize=2048, numbands=1, alpha=0;

		if(numbands>fftsize) {

			numbands = fftsize;
		};

		^this.multiNew('control', fft, fftsize, numbands, alpha);
	}

	init { arg ... theInputs;
		inputs = theInputs;

		^this.initOutputs(theInputs[2], rate);
	}
}