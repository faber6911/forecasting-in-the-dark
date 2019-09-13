<p align="center">
<div align="center">
  <a href="http://datascience.disco.unimib.it/it/"><img src ="https://raw.githubusercontent.com/malborroni/Foundations_of_Computer-Science/master/images/DSunimib.png" width = "100%"></a>
  <h1 align="center">Forecasting in the dark</h1>
  <h6 align="center">a Data Science LAB Project</h6>
<p align="center">
<a href="#overview">Overview &nbsp;</a> |
<a href="#references">&nbsp; References &nbsp;</a> |
<a href="#data">&nbsp; Data &nbsp;</a> |
<a href="#presentation">&nbsp; Presentation &nbsp;</a> |
<a href="#aboutus">&nbsp; About us &nbsp;</a>
  </p>
</div>

<a name="overview"></a>
## &#9741; &nbsp; Overview
<p> <strong>GME</strong> operates in power, gas and environmental markets. It is the exchange place for electricity and natural gas spot trading in Italy.
In the power market platform, producers and purchasers sell and buy wholesale electricity. There is an auction for every hour of the day. 
</p>
Forecasting this <strong>supply function</strong> could be interesting for every energy producer.
<div align = "center">
  <br>
  <br>
  <p align = "center">
  <img align = "center" width = "30%" src = "https://gitlab.com/DBertazioli/enercibiddding/raw/master/img/supply_surface.png">
  </p>
  <p align="center">
  <em>example of supply surface obtained plotting more supply functions all together</em>
  </p>
</div>
<br>
Our objective is to forecast each time series in order to obtain an estimated supply function for future time (1-hour, 24-hours, 168-hours).
<div align = "center">
  <br>
  <br>
  <p align="center">
  <img src = "https://gitlab.com/DBertazioli/enercibiddding/raw/master/img/supply_forecasting.png" width = "35%" align="center">
  </p>
  <p align="center">
  <em>example of time series forecasting</em>
  </p>
</div>
<br>

To do this <strong>Reduced Rank Regression</strong> (RRR) and <strong>Long short-term memory neural nets</strong> (LSTM) are used, in order to compare a pure statistical method with machine learning technique.

Dataset exploration, data preprocessing and LSTM models are executed using Python and all the code is available in the <a href="https://github.com/faber6911/forecasting-in-the-dark/tree/master/PyScripts">PyScripts</a> folder while RRR models are made using R and are available in the <a href="https://github.com/faber6911/forecasting-in-the-dark/tree/master/Rscripts">Rscripts</a> folder.


<a name="references"></a>
## &#9741; &nbsp; References
All the references for this project are available in the <a href="https://github.com/faber6911/forecasting-in-the-dark/tree/master/ref">ref</a> folder.

<a name="data"></a>
## &#9741; &nbsp; Data
All the data used for this project are available in the <a href="https://gitlab.com/DBertazioli/enercibiddding/tree/master/data">data</a> folder (on gitlab) and <a href="https://www.mercatoelettrico.org/it/">here</a>.

<a name="presentation"></a>
## &#9741; &nbsp; Presentation
Our slides presentation is available in the <a href="https://github.com/faber6911/forecasting-in-the-dark/tree/master/slides">slides</a> folder. 

<a name="aboutus"></a>
## &#9741; &nbsp; About us

&#8860; &nbsp; **Alessandro Borroni** 

- *Cosa studio*: Studente Magistrale di Data Science presso l'Università degli Studi di Milano-Bicocca;
- *Studi precedenti*: Laurea triennale in Economia e Commercio presso l'Università degli Studi di Milano-Bicocca.
<br>

<p align = "center">
  <a href = "https://www.linkedin.com/in/alessandro-borroni-212947186/"><img src="https://raw.githubusercontent.com/DBertazioli/Interact/master/img/iconfinder_Popular_Social_Media-22_2329259.png" width = "3%"></a>
  <a href = "https://github.com/malborroni/"><img src="https://raw.githubusercontent.com/malborroni/Foundations_of_Computer-Science/master/images/GitHub.png" width = "3%"></a>
</p>

&#8860; &nbsp; **Dario Bertazioli**

- *Cosa studio*: Studente Magistrale di Data Science presso l'Università degli Studi di Milano-Bicocca;
- *Studi precedenti*: Laurea triennale in Fisica presso l'Università degli Studi di Milano.
<br>

<p align = "center">
  <a href = "https://www.linkedin.com/in/dario-bertazioli-961ab4180/"><img src="https://raw.githubusercontent.com/DBertazioli/Interact/master/img/iconfinder_Popular_Social_Media-22_2329259.png" width = "3%"></a>
  <a href = "https://github.com/DBertazioli/"><img src="https://raw.githubusercontent.com/malborroni/Foundations_of_Computer-Science/master/images/GitHub.png" width = "3%"></a>
</p>

&#8860; &nbsp; **Fabrizio D'Intinosante**

- *Cosa studio*: Studente Magistrale di Data Science presso l'Università degli Studi di Milano-Bicocca;
- *Studi precedenti*: Laurea triennale in Economia e Statistica per le organizzazioni presso l'Università degli Studi di Torino.
<br>

<p align = "center">
  <a href = "https://www.linkedin.com/in/fabrizio-d-intinosante-125042180/"><img src="https://raw.githubusercontent.com/DBertazioli/Interact/master/img/iconfinder_Popular_Social_Media-22_2329259.png" width = "3%"></a>
  <a href = "https://github.com/faber6911/"><img src="https://raw.githubusercontent.com/malborroni/Foundations_of_Computer-Science/master/images/GitHub.png" width = "3%"></a>
</p>

&#8860; &nbsp; **Massimiliano Perletti**

- *Cosa studio*: Studente Magistrale di Data Science presso l'Università degli Studi di Milano-Bicocca;
- *Studi precedenti*: Laurea triennale in Ingegneria dei materiali e delle nano-tecnologie presso il Politecnico di Milano.
<br>

<p align = "center">
  <a href = "https://www.linkedin.com/in/massimilianoperletti/"><img src="https://raw.githubusercontent.com/DBertazioli/Interact/master/img/iconfinder_Popular_Social_Media-22_2329259.png" width = "3%"></a>
  <a href = "https://github.com/maxi93/"><img src="https://raw.githubusercontent.com/malborroni/Foundations_of_Computer-Science/master/images/GitHub.png" width = "3%"></a>
</p>
