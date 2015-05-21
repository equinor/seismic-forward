// $Id: commonheaders.cpp 1279 2014-06-03 14:47:05Z hauge $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// ï  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// ï  Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "commonheaders.hpp"

#include <cstdio>
#include <cstring>
#include "../iotools/fileio.hpp"

using namespace NRLib;
/// Table with EBCDIC-values for the different ASCII values.
/// Values not found in EBCDIC are set to 0x00.
/// Taken from SEG Y rev 1 standard.
const unsigned char transtab[] =
// 0     1     2     3     4     5     6     7
{0x00, 0x01, 0x02, 0x03, 0x37, 0x2d, 0x2e, 0x2f,
 0x16, 0x05, 0x25, 0x0b, 0x0c, 0x0d, 0x00, 0x00,
 0x10, 0x11, 0x12, 0x13, 0x3c, 0x3d, 0x32, 0x26,
 0x18, 0x00, 0x3f, 0x27, 0x1c, 0x1d, 0x1e, 0x1f,
 0x40, 0x5a, 0x7f, 0x7b, 0x5b, 0x6c, 0x50, 0x7d,
 0x4d, 0x5d, 0x5c, 0x4e, 0x6b, 0x60, 0x4b, 0x61,
 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7,
 0xf8, 0xf9, 0x7a, 0x5e, 0x4c, 0x7e, 0x6e, 0x6f,
 0x7c, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7,
 0xc8, 0xc9, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6,
 0xd7, 0xd8, 0xd9, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6,
 0xe7, 0xe8, 0xe9, 0x00, 0xe0, 0x00, 0x00, 0x6d,
 0x79, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
 0x88, 0x89, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96,
 0x97, 0x98, 0x99, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6,
 0xa7, 0xa8, 0xa9, 0xc0, 0x4f, 0xd0, 0xa1, 0x07
 };

/// Translates an ASCII character to EBCDIC.
static unsigned char ascii2ebcdic(unsigned char c)
{
  if (c > 128)
    return 0;
  return transtab[c];
}

static void emptyHeader(std::string &ebcdicH)
{
  ebcdicH = "";
   ebcdicH+=  "√@Ò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+=  "√@Ú@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@Û@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@Ù@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@ı@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@ˆ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@˜@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@¯@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@˘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√Ò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√ÒÒ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√ÒÚ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√ÒÛ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÒÙ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Òı@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Òˆ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ò˜@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ò¯@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ò˘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ú@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÚÒ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÚÚ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÚÛ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÚÙ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Úı@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Úˆ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ú˜@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ú¯@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ú˘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Û@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÛÒ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÛÚ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÛÛ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÛÙ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ûı@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ûˆ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Û˜@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Û¯@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Û˘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ù@≈’ƒ@≈¬√ƒ…√@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+='@';
}

static void oldHeader(std::string &ebcdicH)
{
  ebcdicH = "";
  ebcdicH+= "√@Ò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@Ú@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@Û@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@Ù@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@ı@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@ˆ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@˜@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@¯@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√@˘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√Ò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√ÒÒ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√ÒÚ@@@@@@@@@„»…‚@‚≈«Ë@∆…”≈@Ê¡‚@ÊŸ…„„≈’@¬Ë@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√ÒÛ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√ÒÙ@@@@@@@@@„»≈@√Ÿ¡Â¡@◊Ÿ÷«Ÿ¡‘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+= "√Òı@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Òˆ@@@@@@@@@ƒ≈Â≈”÷◊≈ƒ@¬Ë@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ò˜@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ò¯@@@@@@@@@Ÿ¡«’¡Ÿ@»¡‰«≈@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ò˘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ú@@@@@@@@@¡’’≈@Ÿ¡’ƒ…@‚ËÂ≈Ÿ‚Â≈≈’@@@@@@@@@@@¡’ƒ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÚÒ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÚÚ@@@@@@@@@÷ƒƒ@“÷”¬—÷Ÿ’‚≈’@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÚÛ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÚÙ@@@@@@@@@¡„@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Úı@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Úˆ@@@@@@@@@„»≈@’÷ŸÊ≈«…¡’@√÷‘◊‰„…’«@√≈’„≈Ÿ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ú˜@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ú¯@@@@@@@@@ÚÛ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ú˘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Û@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÛÒ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÛÚ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÛÛ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√ÛÙ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ûı@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ûˆ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Û˜@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Û¯@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Û˘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="√Ù@≈’ƒ@≈¬√ƒ…√@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  ebcdicH+="@";
}

TextualHeader::TextualHeader()
{
  emptyHeader(buffer_);
  //buffer_[3199] = '@';
}

TextualHeader TextualHeader::standardHeader()
{
  TextualHeader header;
  oldHeader(header.buffer_);
  return header;
}

int TextualHeader::SetLine(int lineNo, const std::string& text)
{
  int retCode = 0;

  // Don't allow modification of last line.
  if (lineNo < 0 || lineNo > 39) {
    return -1;
  }

  size_t length = text.size();
  if (length > 75) {
    length = 75;
    retCode = 1;
  }

  size_t pos = lineNo*80 + 5;
  for (size_t i = 0; i < length; ++i) {
    buffer_[pos + i] = ascii2ebcdic(text[i]);
  }
  for (size_t i = length; i < 75; ++i) {
    buffer_[pos + i] = ascii2ebcdic(' ');
  }

  return retCode;
}


void TextualHeader::Write(std::ostream& file) const
{
  //fwrite(buffer_, 1, 3200, file);
 // char buffer[3200];
//  int i;
 // for (i=0;i<3200;i++)
 //   buffer[i] = buffer_[i];
  for (size_t i=0;i<3200;i++)
    file << buffer_[i] ;
  //file.write(buffer_,3200);
}


BinaryHeader::BinaryHeader(std::istream& file)
{
  Update(file);
}


BinaryHeader::BinaryHeader()
{
}


void
BinaryHeader::Update(std::istream& file)
{
  ReadBinaryInt(file);
  lino_ = ReadBinaryInt(file);
  ReadBinaryInt(file);
  ReadBinaryShort(file);
  ReadBinaryShort(file);
  hdt_ = ReadBinaryShort(file); ///=dz
  ReadBinaryShort(file);
  hns_ = ReadBinaryShort(file); ///=nz
  ReadBinaryShort(file);
  format_ = ReadBinaryShort(file);

  char buffer[374];
  file.read(buffer,374);
}

void BinaryHeader::Write(std::ostream& file, double dz, size_t nz, short n_sam_per_ens) const
{
  int dummy = 9999;
  WriteBinaryInt(file,dummy); //1-4
  dummy = 1;
  WriteBinaryInt(file,dummy); //5-8
  WriteBinaryInt(file,dummy); //9-12

  WriteBinaryShort(file, n_sam_per_ens); //13-14
  WriteBinaryShort(file, n_sam_per_ens); //15-16
  WriteBinaryShort(file, static_cast<short>(1000*dz)); // Hdt 17-18
  short dummy2 = 0;
  WriteBinaryShort(file, dummy2); //19-20
  WriteBinaryShort(file, static_cast<short>(nz)); ///=nz Hns 21-22
  WriteBinaryShort(file, dummy2); //23-24
  dummy2 = 1;
  WriteBinaryShort(file, dummy2); // format 25-26
  WriteBinaryShort(file, n_sam_per_ens); // 27-28
  dummy2 = 4;
  WriteBinaryShort(file, dummy2); //29-30
  dummy2 = 0;
  int i;
  for (i=0;i<12;i++)
    WriteBinaryShort(file, dummy2); //31-54
  dummy2 = 1;
  WriteBinaryShort(file, dummy2); //55-56
  char buffer[246];
  memset(buffer, 0, 246);
  file.write(buffer,246); //57-302
  WriteBinaryShort(file, dummy2); //Constant tracelength, 303-304
  file.write(buffer, 96); //305-400
}
