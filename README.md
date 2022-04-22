# A Novel Algorithm for Heart Rate Estimation Based on Synchrosqueezing Transform [1]

To estimate the heart rate (HR) from electrocardiogram (ECG), time-frequency representations such as the short-time Fourier transform (STFT) is often used. As the STFT is constrained by the choice of a specific analysis window used for its definition, we alternatively propose to estimate HR from a synchrosqueezed STFT. More precisely, we build a novel algorithm inspired by non-negative matrix factorization to estimate HR by determining the minimal Wasserstein distance between a synchrosqueezed STFT and columns of a specific dictionary matrix. Throughout numerical simulations carried out on both synthetic and real ECGs, we show in what way to use a synchrosqueezed STFT rather than STFT improves HR estimation.

[1] Nils Laurent,Sylvain Meignen,Julie Fontecave-Jallon,Bertrand Rivet -- A Novel Algorithm for Heart Rate Estimation Based on Synchrosqueezing Transform -- EUSIPCO 2021 - 29th European Signal Processing Conference (EUSIPCO), pp. 1286-1290

## Implementation

It relies on:
* [2] the 'NMF toolbox': https://www.audiolabs-erlangen.de/resources/MIR/NMFtoolbox/
* [3] the synthetic ECG generator : https://physionet.org/files/ecgsyn/1.0.0/Matlab/

[2] Patricio López-Serrano, Christian Dittmar, Yiğitcan Özer, and Meinard Müller NMF Toolbox: Music Processing Applications of Nonnegative Matrix Factorization In Proceedings of the International Conference on Digital Audio Effects (DAFx), 2019.

[3] Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.

## Data

Generate the data using `gen_data_ecgsyn_N.m` and `gen_data_ecgsyn_single.m`
