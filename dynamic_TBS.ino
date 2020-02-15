#include <TimerOne.h>

// This code was used for the project -
// "The effects of individualised intermittent theta burst stimulation in the prefrontal cortex: A TMS‐EEG study"
// for varying theta-burst stimulation parameters

//Input Variables
                             // Warp time for Debug
                         // Increment Time every 10us
int nP=3;                                   // Number of Pulses in a Burst
int nB= 10;                                 // Number of Bursts in a Train
int nT = 20;                                // Number of Trains in a Session
//float sF=1/200;
//float pF=1/50;
//float bF=1/5;
//float tF=8/1;
int timestep=10;
unsigned long sF=5;
unsigned long pF=20;//20; //20;//33  ; gamma
unsigned long  bF=200;//200; //200;//167; theta
unsigned long tF=8000;//10000;
 int warp=1000;  
 int offset=62;
unsigned long onesec=1000000/(timestep*warp);
unsigned long Pwidth = 500;                // trigger duration
unsigned long IPI = (((pF)*onesec));// (-Pwidth)+(((1/50)*onesec)/;               // Inter-Pulse-Interval
unsigned long PDur =(IPI+Pwidth);              // Pulse Duration {IPI + Pwidth}
unsigned long IBI =((bF)*onesec);// ((bF)*onesec)-(2*PDur)/2;//-(2*PDur)+(((1/5)*onesec)/timestep);              // Inter-Burst-Inerval
unsigned long BDur =nP*(PDur)  ;            // Burst Duration {nP*(IBI+PDur)}
unsigned long ITI =(((tF)*onesec));// (((tF)*onesec)/2)+IBI+BDur+IPI;//+9*IBI;// 80000*warp;// Inter-Train-Interval
unsigned long TDur =nB*(BDur+IBI) ;           // Train Duration {nB*(ITI+BDur)}//unsigned long IP1 = 100*warp;               // for biphasic first pulse width

//unsigned long IP1 = 100*warp;               // for biphasic first pulse width
//unsigned long IP2 = 0*warp;                 // for biphasic  second pulse width
// Output Configuration
const int ledPin =  12; //  4;   // the number of the LED pin
const int triggerPin =3;//  4;    // the number of the Trigger pin


// Internal Variables
unsigned long SystemTime;    
unsigned long t1;
unsigned long t2;
unsigned long elapsed;  
unsigned long realtime;              
unsigned long LastLoopTime;          
int state = 11;
int long PulseTrainTimestamps;
int long PrePulseTrainTimestamps;
int StimulatingState = 0;
int            PreStimulusStatus = 0;
int           StimulusStatus = 0;
int           PulseStatus = 0;
int            BurstStatus = 0;
unsigned long StimulusTrainEndTime;
unsigned long NextBurstTransitionTime = 0;
unsigned long NextPulseTransitionTime = 0;
unsigned long nextTime= 0;
int Tcount ;
int Bcount;
int Pcount;
int Ptotal;
int stage=1;
byte CommandByte;
byte MessageByte;
int flip;
void setup() {
  Serial.begin(115200);
  pinMode(ledPin, OUTPUT); // set the digital pin as output:
  pinMode(triggerPin, OUTPUT);
  SystemTime = 0;
  LastLoopTime = SystemTime;
  Timer1.attachInterrupt(handler);
  Timer1.initialize(timestep);               // Calls every 10us ;
  Timer1.start();
}
void loop() {
 // //time goes by....tick tock
//Serial.println(Pcount);
//Serial.println(0);
//Serial.println(stage);
//Serial.println(nextTime);
}

void handler(void) {
  //interupt evey [timestep]microseconds
  SystemTime++;
//realtime =realtime+micros();
 //realtime =0+micros();
// realtime++;
//  if (Serial.available() >= 2) {  // wait for serial communication
//    // read the incoming byte:
//    CommandByte = Serial.read(); // listen byte {213}
//    if (CommandByte == 'c') {
//      CommandByte = Serial.read();              // switch case byte {72 = handshake ; 74= read param}
//       MessageByte = Serial.read();
//      Serial.println(CommandByte);
//      //state =1;
//      StimulatingState = 1;
//      PrePulseTrainTimestamps = SystemTime;
//      state = 1;
//    }
//  }

  //_____________
  if (SystemTime >= (PrePulseTrainTimestamps + nextTime)) {
PrePulseTrainTimestamps =0;
//Serial.println(nextTime);
//Serial.println(SystemTime);
   
//Serial.println(elapsed);

    //Serial.println(state);
//Serial.println(stage*SystemTime);
//Serial.println(Bcount);
//Serial.println(Tcount);
//Serial.println(0);
//Serial.println(SystemTime);


//____________
 //////////////////////////////////

switch (state){
case 0:
nextTime =0;
state=0;
break;

case 1:
nextTime=ITI/2;
//digitalWrite(ledPin, HIGH);
state=4;
break;

case 2:
nextTime =IBI/2;
//digitalWrite(ledPin, LOW);
state=4;
break;

case 3:
nextTime =IPI/2;
state=4;
break;

case 4:
nextTime =Pwidth;
dacWriteON();
//Serial.println(realtime);
realtime=0;
PrePulseTrainTimestamps=0;
//Serial.println(Pwidth);
state=44;
break;

case 44:
dacWriteOFF();
Ptotal++;
///////////////////////////////
if (Pcount<nP){
       Pcount=Pcount;
       Bcount=Bcount;
       Pcount++;
        nextTime =IPI-Pwidth;//((IPI-Pwidth)/2)+(Pwidth/2)-155;
        //nextTime=1000;// LOCKED 0.2SEC interval!!!!
        //Serial.println(nextTime);
        //Serial.println(Ptotal);
        state=4;
        PrePulseTrainTimestamps =0;//offset;//-64;
 }/////////////////////////
if (Pcount>=nP){
      Pcount=0;
      Bcount=Bcount;
      Bcount++;
      state=4;
      nextTime =IBI-(IPI+PDur);//IBI/2 +BDur+500+IPI;//IBI-(IPI+PDur);////-(IPI+500+500);//+(IPI-(nP*Pwidth));
      // nextTime=10250;// LOCKED 0.2SEC interval!!!!
      PrePulseTrainTimestamps =0;//-offset;
//     Serial.println(IPI);
//      Serial.println(IBI);
      //nextTime=0;
// }else{
//    Pcount++;
////    state=state;
 }/////////////////////////
if (Bcount>nB-1){
            Pcount=0;
            Bcount=0;
            Tcount++;
            state=4;
            nextTime=ITI;// 522500;// LOCKED 8SEC interval!!!!
            //nextTime =ITI+(nB*(IBI+((nP*IPI)-(nP*Pwidth))));//ITI+(((BDur+(nP*PDur))*nB)/2);//*2;//+(TDur/2)+IBI;//+((nP*64)*nB);//+BDur+PDur+Pwidth;//IBI+IPI+Pwidth;//+(3*(9*(IBI-(2*(PDur)))));
            PrePulseTrainTimestamps =0;//(((nP)*offset)*nB)*2;//((nP*offset)*nB);//;nB(*64;//(nB-1)*(nP*64);//-((9*(BDur)));//-((64)*10000); //((64*nB)+(nP*64));//50;
            //Serial.println(nextTime);
            //Serial.println(Tcount);
            //nextTime=0;
//          }else{
//            state=state;
}///////////////////////////
if(Tcount>=nT){
              Pcount=0;
              Bcount=0;
              Tcount=0;
              state=0;
              state=0;
//              }else{
//              state=state;
}////////////////////////////
//                }
//          }
//    }

break;
case 33:
break;
case 22:
nextTime=10000*onesec/2; //
digitalWrite(ledPin, LOW);
state=1;
break;
case 11:
nextTime=10000*onesec; //
state=22;
digitalWrite(ledPin, HIGH);
break;
}
//////////////////////////////////
  //_____________

SystemTime = 0;

//elapsed=t2-t1;
//////////////////////////////////

}
 
}//END handler


// FUNCTIONS //
void dacWriteON() {
      digitalWrite(ledPin, HIGH);
    digitalWrite(triggerPin, HIGH);
   //Serial.println(realtime);
}

void dacWriteOFF() {
   
    digitalWrite(triggerPin, LOW);
    digitalWrite(ledPin, LOW);
   
    //realtime=0;
}

//void ITdelay() {
//   Timer1.stop();
//    delay(8000);
//    Timer1.start();
//    //realtime=0;
//}


//void dacWrite(int input) {
//  if (input == 1) {
//      digitalWrite(ledPin, HIGH);
//    digitalWrite(triggerPin, HIGH);
//  }else{
//    digitalWrite(triggerPin, LOW);
//    digitalWrite(ledPin, LOW);
//    }
//  }




//END FUNCTIONS //

//
//if (StimulatingState == 0) {
//  SystemTime = 0;



//
// if (SystemTime >= (PrePulseTrainTimestamps + IPI)) {
//break;
//StimulusTrainEndTime = SystemTime + TDur;
//NextBurstTransitionTime = SystemTime+BDur;
