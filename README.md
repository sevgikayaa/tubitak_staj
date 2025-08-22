# tubitak_staj

11 Ağustos 2025 - 19 Eylül 2025 tarihleri arasında Gebze TÜBİTAK Temel Bilimler Araştırma Enstitüsü'nde yaptığım gönüllü staj süresi boyunca yaptıklarım.

gaia_hrdiagram.ipynb adlı dosyada yaptığım şey astroquery ile gaia'dan veri çekerek interaktif bir HR diyagramı oluşturmak. HR diyagramı oluştuğunda yıldızların üstüne geldiğinizde yıldızın ra, dec, magnitude, color index ve source id bilgilerini gösteriyor. Renklendirme olarak daha parlak yıldızlar sarı gözükürken daha sönük yıldızlar mor gözükmekte (görsellik açısından). Toplam 36000 yıldız verisi çekmiş oldu, bu yıldızlar belli kriterlere uygun yıldızlar bu kriterler araştırmamıza göre değiştirilebilir.

galaktik_duzlem.ipynb adlı dosyada güneş merkezde olacak şekilde, güneş çevresinde bulunan belli kriterlerdeki yıldızların galaktik koordinatlarına göre dağılımı ve yoğunluğu 3D olacak şekilde verilmiştir. Gaia veritabanından yakın yıldızlar çekilmiş olup, yıldızların RA/Dec ve paralaks bilgileri ile galaktik konumları (x,y,z) hesaplanmıştır. Radius parametresi ile uzaklığı interaktif değiştirerek yıldız yoğunluğu gözlenebilmektedir. toplam 100000 yıldız verisi çekilmiştir.

exoplanet_transit.ipynb adlı dosyada Nasa exoplanet archive üzerinden ötegezegenlerin verisi, TESS üzerinden ise bu ötegezegenlerin yıldızlarının verileri (ışık eğrileri) çekilerek bir ötegezegen ışık eğrisi oluşturulmuştur. Test amaçlı yalnızca 1 adet ötegezegenin eğrisi bulunmaktadır. 1 adet verinin bile işlenmesi yaklaşık 1 saat sürmüştür. Limit değeri None yazılarak keşfedilen bütün ötegezegenlerin eğrisi çizdirilebilir.

lightkurve_data_analysis adlı doyada lightkurve kütüphanesi ile MAST portal üzerinden seçilen ötegezegen analizi ve TESS üzerinden gelen ötegezegen verilerini (pikselleri) ışık eğrisine çevirerek analizi yapılmıştır.

KIC-11026764.ipynb ve KIC-9832227_KIC-11145123.ipynb --> ötegezegeni bulunmayan veya transit göstermeyen yıldızların ışık eğrisi analizi sonucu elde edilen grafikler
Kepler-42-d.ipynb --> Kepler-42 yıldızının Kepler-42 d ötegezegeninin ışık eğrisi analizi
