UInt_t SetSeedOwn(){
    UInt_t  fSeed;
    UInt_t  fSeed1;
    UInt_t  fSeed2;
      TUUID u; 
      UChar_t uuid[16];
      u.GetUUID(uuid);
      fSeed  =  UInt_t(uuid[3])*16777216 + UInt_t(uuid[2])*65536 + UInt_t(uuid[1])*256 + UInt_t(uuid[0]);
      fSeed1  =  UInt_t(uuid[7])*16777216 + UInt_t(uuid[6])*65536 + UInt_t(uuid[5])*256 + UInt_t(uuid[4]);
      fSeed2  =  UInt_t(uuid[11])*16777216 + UInt_t(uuid[10])*65536 + UInt_t(uuid[9])*256 + UInt_t(uuid[8]);
      // use also the other bytes
      UInt_t seed3 = UInt_t(uuid[15])*16777216 + UInt_t(uuid[14])*65536 + UInt_t(uuid[13])*256 + int(uuid[12]);
      fSeed2 += seed3; 



      //    TRandom r3(0);
      // fSeed   = static_cast<UInt_t> (4294967296.*r3.Rndm());
      // fSeed1  = static_cast<UInt_t> (4294967296.*r3.Rndm());
      // fSeed2  = static_cast<UInt_t> (4294967296.*r3.Rndm());

      if (fSeed < 2)   fSeed  += 2UL;
      if (fSeed1 < 8)  fSeed1 += 8UL;
      if (fSeed2 < 16) fSeed2 += 16UL;

     UInt_t iy = fSeed ^ fSeed1 ^ fSeed2;
     if(iy)       return iy;
     else return 0;
}
