// Minimal STORE-method ZIP writer (no compression). Sufficient for bundling
// agent-emitted output files (GenBank, CSV, JSON, TXT) on the client without
// pulling in JSZip. Format reference: PKZIP APPNOTE 4.5 sections 4.3.7 (local
// file header), 4.3.12 (central directory), 4.3.16 (end of central directory).

export type ZipEntry = { name: string; bytes: Uint8Array };

const CRC32_TABLE: Uint32Array = (() => {
  const t = new Uint32Array(256);
  for (let i = 0; i < 256; i++) {
    let c = i;
    for (let k = 0; k < 8; k++) {
      c = (c & 1) ? (0xedb88320 ^ (c >>> 1)) : (c >>> 1);
    }
    t[i] = c >>> 0;
  }
  return t;
})();

function crc32(bytes: Uint8Array): number {
  let c = 0xffffffff;
  for (let i = 0; i < bytes.length; i++) {
    c = CRC32_TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
  }
  return (c ^ 0xffffffff) >>> 0;
}

function utf8(s: string): Uint8Array {
  return new TextEncoder().encode(s);
}

// MS-DOS timestamp (seconds resolution / 2, since 1980). Used in both local and
// central headers.
function dosTime(d: Date): { time: number; date: number } {
  const time = ((d.getHours() & 0x1f) << 11) |
               ((d.getMinutes() & 0x3f) << 5) |
               (((d.getSeconds() / 2) | 0) & 0x1f);
  const date = (((d.getFullYear() - 1980) & 0x7f) << 9) |
               (((d.getMonth() + 1) & 0x0f) << 5) |
               (d.getDate() & 0x1f);
  return { time, date };
}

export function writeZip(entries: ZipEntry[], now: Date = new Date()): Blob {
  const { time, date } = dosTime(now);
  const parts: Uint8Array[] = [];
  const central: Uint8Array[] = [];
  let offset = 0;

  for (const entry of entries) {
    const nameBytes = utf8(entry.name);
    const data = entry.bytes;
    const crc = crc32(data);

    // Local file header (30 bytes + name + extra)
    const lh = new Uint8Array(30 + nameBytes.length);
    const lhv = new DataView(lh.buffer);
    lhv.setUint32(0, 0x04034b50, true);     // local file header signature
    lhv.setUint16(4, 20, true);              // version needed (2.0)
    lhv.setUint16(6, 0x0800, true);          // general purpose flag (UTF-8)
    lhv.setUint16(8, 0, true);               // compression method (STORE)
    lhv.setUint16(10, time, true);           // mod file time
    lhv.setUint16(12, date, true);           // mod file date
    lhv.setUint32(14, crc, true);            // CRC-32
    lhv.setUint32(18, data.length, true);    // compressed size
    lhv.setUint32(22, data.length, true);    // uncompressed size
    lhv.setUint16(26, nameBytes.length, true); // file name length
    lhv.setUint16(28, 0, true);              // extra field length
    lh.set(nameBytes, 30);

    parts.push(lh);
    parts.push(data);

    // Central directory file header (46 bytes + name + extra + comment)
    const ch = new Uint8Array(46 + nameBytes.length);
    const chv = new DataView(ch.buffer);
    chv.setUint32(0, 0x02014b50, true);      // central file header signature
    chv.setUint16(4, 0x0314, true);          // version made by (UNIX, 2.0)
    chv.setUint16(6, 20, true);              // version needed (2.0)
    chv.setUint16(8, 0x0800, true);          // general purpose flag (UTF-8)
    chv.setUint16(10, 0, true);              // compression method
    chv.setUint16(12, time, true);
    chv.setUint16(14, date, true);
    chv.setUint32(16, crc, true);
    chv.setUint32(20, data.length, true);
    chv.setUint32(24, data.length, true);
    chv.setUint16(28, nameBytes.length, true);
    chv.setUint16(30, 0, true);              // extra field length
    chv.setUint16(32, 0, true);              // file comment length
    chv.setUint16(34, 0, true);              // disk number start
    chv.setUint16(36, 0, true);              // internal file attributes
    chv.setUint32(38, 0, true);              // external file attributes
    chv.setUint32(42, offset, true);         // relative offset of local header
    ch.set(nameBytes, 46);

    central.push(ch);
    offset += lh.length + data.length;
  }

  const centralStart = offset;
  let centralSize = 0;
  for (const c of central) centralSize += c.length;

  // End of central directory record
  const eocd = new Uint8Array(22);
  const ev = new DataView(eocd.buffer);
  ev.setUint32(0, 0x06054b50, true);          // EOCD signature
  ev.setUint16(4, 0, true);                    // disk number
  ev.setUint16(6, 0, true);                    // disk where central dir starts
  ev.setUint16(8, entries.length, true);       // central dir records on this disk
  ev.setUint16(10, entries.length, true);      // total central dir records
  ev.setUint32(12, centralSize, true);         // size of central directory
  ev.setUint32(16, centralStart, true);        // offset of central dir start
  ev.setUint16(20, 0, true);                   // comment length

  const blobParts: BlobPart[] = [...parts, ...central, eocd] as BlobPart[];
  return new Blob(blobParts, { type: "application/zip" });
}
