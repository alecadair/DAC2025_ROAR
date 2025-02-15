v {xschem version=3.4.5 file_version=1.2
}
G {}
K {}
V {}
S {}
E {}
N 80 -360 300 -360 {
lab=avdd_1v8}
N 280 -280 300 -280 {
lab=GND}
N 500 -320 550 -320 {
lab=out}
N 550 -260 550 -240 {
lab=GND}
N 400 -270 400 -250 {
lab=iref}
N 170 -340 300 -340 {
lab=inn}
N 80 -340 80 -280 {
lab=inn}
N 240 -300 240 -240 {
lab=inp}
N 240 -300 300 -300 {
lab=inp}
N -10 -360 80 -360 {
lab=avdd_1v8}
N 80 -340 170 -340 {
lab=inn}
C {devices/vsource.sym} -10 -330 0 0 {name=V1 value=1.8 savecurrent=false}
C {devices/gnd.sym} -10 -300 0 0 {name=l1 lab=GND}
C {devices/gnd.sym} 280 -280 0 0 {name=l2 lab=GND}
C {devices/capa.sym} 550 -290 0 0 {name=C1
m=1
value=50f
footprint=1206
device="ceramic capacitor"}
C {devices/gnd.sym} 550 -240 0 0 {name=l3 lab=GND}
C {devices/gnd.sym} 400 -190 0 0 {name=l4 lab=GND}
C {devices/vsource.sym} 80 -250 0 0 {name=V3
value="SIN(0.9 0.1 1Meg)"}
C {devices/gnd.sym} 80 -220 0 0 {name=l5 lab=GND}
C {devices/vsource.sym} 240 -210 0 0 {name=V4 value=0.9 savecurrent=false}
C {devices/gnd.sym} 240 -180 0 0 {name=l6 lab=GND}
C {devices/lab_wire.sym} 540 -320 0 0 {name=p2 sig_type=std_logic lab=out}
C {devices/lab_pin.sym} 270 -340 0 0 {name=p3 sig_type=std_logic lab=inn}
C {devices/lab_pin.sym} 280 -300 0 0 {name=p4 sig_type=std_logic lab=inp}
C {devices/code.sym} 560 -610 0 0 {name=TT_MODELS
only_toplevel=true
format="tcleval( @value )"
value="
** opencircuitdesign pdks isntsall
.lib /opt/pdk/open_pdks/sky130/sky130A/libs.tech/ngspice/sky130.lib.spice tt

"
spice_ignore=false}
C {devices/code_shown.sym} 720 -670 0 0 {name="AC Analysis" only_toplevel=false value="

.temp 27
.ac dec 100 10 100G
.print ac V(out)
.measure ac unity_freq when vdb(out)=0
//.measure ac gain_max MAX vdb(out)
.measure ac pole1_freq when vp(out)=-45 cross=1
.measure ac pole2_freq when vp(out)=-90 cross=1
//.measure ac bandwidth when vdb(out)= '(gain_max - 3)' cross=1
.control
	//run
	//save all
	//ac dec 10 1 1e11
	//settype decibel out
	//setcurplottitle=Loopgain
	//let frequency ac1.frequency
.endc

"}
C {cm_ota.sym} 300 -260 0 0 {name=X1 w1_2=0.950 w3_4=0.630 w5_6=3.78 w7_8=2.52}
C {devices/code.sym} 630 -290 0 0 {name=Transient only_toplevel=false value="

.temp 27
.option savecurrents
.save all
.control
	op
	write tb_cm_ota_working_v2_extracted.raw
	tran 0.1ns 10us
	plot v(inn) v(out)
.endc
"}
C {devices/lab_pin.sym} 400 -260 2 0 {name=p5 sig_type=std_logic lab=iref
}
C {devices/lab_pin.sym} -10 -360 0 0 {name=p6 sig_type=std_logic lab=avdd_1v8}
C {devices/isource.sym} 400 -220 2 0 {name=I1 value=24u}
