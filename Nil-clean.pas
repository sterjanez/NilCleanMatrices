program nilclean;
uses
  crt;
const
  nmax=20;
  maxpower=16;
function blocno_zgornje_trikotni_idempotenti(ime1,ime2,ime3:string):longint;
  label 1,2;
  var
    dat1,dat2,dat3:text;
    idempotenti:longint;
    i,j,k,m,n:byte;
    e:array[1..nmax,1..nmax] of boolean;
    vsota:boolean;
    s:string;
  begin
    assign(dat1,ime1);
    assign(dat2,ime2);
    assign(dat3,ime3);
    rewrite(dat3);
    reset(dat1);
    idempotenti:=0;
    while not eof(dat1) do
    begin
      readln(dat1,s);
      m:=round(sqrt(length(s)));
      for i:=1 to m do
        for j:=1 to m do
        begin
          e[i,j]:=s[1]='1';
          delete(s,1,1);
        end;
      reset(dat2);
      while not eof(dat2) do
      begin
        readln(dat2,s);
        n:=round(sqrt(length(s)));
        for i:=m+1 to m+n do
          for j:=m+1 to m+n do
          begin
            e[i,j]:=s[1]='1';
            delete(s,1,1);
          end;
        for i:=1 to m do
          for j:=m+1 to m+n do
            e[i,j]:=false;
        for i:=m+1 to m+n do
          for j:=1 to m do
            e[i,j]:=false;
        2:
        for i:=1 to m do
          for j:=m+1 to m+n do
          begin
            vsota:=false;
            for k:=1 to m+n do
              vsota:=vsota xor (e[i,k] and e[k,j]);
            if vsota<>e[i,j] then
              goto 1;
          end;
        s:='';
        for i:=1 to m+n do
          for j:=1 to m+n do
            if e[i,j] then
              s:=s+'1'
            else
              s:=s+'0';
        writeln(dat3,s);
        inc(idempotenti);
        1:
        for i:=1 to m do
          for j:=m+1 to m+n do
            if e[i,j] then
              e[i,j]:=false
            else
            begin
              e[i,j]:=true;
              goto 2;
            end;
      end;
      close(dat2);
    end;
    close(dat1);
    close(dat3);
    blocno_zgornje_trikotni_idempotenti:=idempotenti;
  end;
function izlusci_idempotente(s,ime1,ime2:string):longint;
  label 1;
  var
    n,i,j,k:byte;
    a,q,r:array[1..nmax,1..nmax] of boolean;
    vsota:boolean;
    dat1,dat2:text;
    z:string;
    idempotenti:longint;
  begin
    n:=round(sqrt(length(s)));
    for i:=1 to n do
      for j:=1 to n do
      begin
        a[i,j]:=s[1]='1';
        delete(s,1,1);
      end;
    assign(dat1,ime1);
    assign(dat2,ime2);
    reset(dat1);
    rewrite(dat2);
    idempotenti:=0;
    while not eof(dat1) do
    begin
      readln(dat1,s);
      z:=s;
      for i:=1 to n do
        for j:=1 to n do
        begin
          q[i,j]:=a[i,j] xor (s[1]='1');
          delete(s,1,1);
        end;
      for i:=1 to n do
        for j:=1 to n do
        begin
          vsota:=false;
          for k:=1 to n do
            vsota:=vsota xor (q[i,k] and q[k,j]);
          r[i,j]:=vsota;
        end;
      for i:=1 to n do
        for j:=1 to n do
        begin
          vsota:=false;
          for k:=1 to n do
            vsota:=vsota xor (q[i,k] and r[k,j]);
          if vsota then
            goto 1;
        end;
      writeln(dat2,z);
      inc(idempotenti);
      1:
    end;
    close(dat1);
    close(dat2);
    izlusci_idempotente:=idempotenti;
  end;
function zapisi_bzt_idempotente(s,ime:string):longint;
  label 1,2;
  var
    z:string;
    n,r,i,j,k,n0:byte;
    ni:array[1..nmax] of byte;
    e,e2minuse,stevec:array[1..nmax,1..nmax] of boolean;
    idempotenti:longint;
    dat:text;
  begin
    s:=s+',';
    z:=s;
    n:=0;
    r:=0;
    while z<>'' do
    begin
      r:=r+1;
      i:=pos(',',z);
      ni[r]:=round(sqrt(i-1));
      n:=n+ni[r];
      delete(z,1,i);
    end;
    for i:=1 to n do
      for j:=1 to n do
      begin
        e[i,j]:=false;
        e2minuse[i,j]:=false;
        stevec[i,j]:=false;
      end;
    n0:=0;
    for k:=1 to r do
    begin
      for i:=n0+1 to n0+ni[k] do
        for j:=n0+1 to n0+ni[k] do
        begin
          e[i,j]:=s[1]='1';
          delete(s,1,1);
        end;
      delete(s,1,1);
      n0:=n0+ni[k];
    end;
    assign(dat,ime);
    rewrite(dat);
    idempotenti:=0;
    1:
    for i:=1 to n do
      for j:=1 to n do
        if e2minuse[i,j] then
          goto 2;
    s:='';
    for i:=1 to n do
      for j:=1 to n do
        if e[i,j] then
          s:=s+'1'
        else
          s:=s+'0';
    writeln(dat,s);
    inc(idempotenti);
    2:
    n0:=0;
    for k:=1 to r do
    begin
      for i:=n0+1 to n0+ni[k] do
        for j:=n0+ni[k]+1 to n do
          if stevec[i,j] then
            stevec[i,j]:=false
          else
          begin
            stevec[i,j]:=true;
            for k:=1 to n do
            begin
              if e[k,i] then
                e2minuse[k,j]:=not e2minuse[k,j];
              if e[j,k] then
                e2minuse[i,k]:=not e2minuse[i,k];
            end;
            if i<>j then
            e2minuse[i,j]:=not e2minuse[i,j];
            e[i,j]:=not e[i,j];
            goto 1;
          end;
      n0:=n0+ni[k];
    end;
    close(dat);
    zapisi_bzt_idempotente:=idempotenti;
  end;
function zapisi_idempotente(n:byte;ime_datoteke:string):longint;
  label 1,2,3;
  var
    i,j,k,nplus1:byte;
    e,e2minuse,stevec:array[1..nmax+1,1..nmax+1] of boolean;
    s:string;
    idempotenti:longint;
    dat:text;
  begin
    nplus1:=n+1;
    for i:=1 to nplus1 do
      for j:=1 to nplus1 do
      begin
        e[i,j]:=false;
        e2minuse[i,j]:=(i=nplus1) or (j=nplus1);
        stevec[i,j]:=false;
      end;
    assign(dat,ime_datoteke);
    rewrite(dat);
    idempotenti:=0;
    1:
    i:=1;
    j:=1;
    2:
    if not e2minuse[i,j] then
    begin
      if j=n then
      begin
        j:=1;
        inc(i);
      end
      else
        inc(j);
      goto 2;
    end;
    if i=nplus1 then
    begin
      s:='';
      for i:=1 to n do
        for j:=1 to n do
          if e[i,j] then
            s:=s+'1'
          else
            s:=s+'0';
      writeln(dat,s);
      inc(idempotenti);
    end;
    i:=1;
    j:=1;
    3:
    if stevec[i,j] then
    begin
      stevec[i,j]:=false;
      if j=n then
      begin
        j:=1;
        inc(i);
      end
      else
        inc(j);
      goto 3;
    end;
    if i<>nplus1 then
    begin
      stevec[i,j]:=true;
      for k:=1 to n do
      begin
        if e[k,i] then
          e2minuse[k,j]:=not e2minuse[k,j];
        if e[j,k] then
          e2minuse[i,k]:=not e2minuse[i,k];
      end;
      if i<>j then
        e2minuse[i,j]:=not e2minuse[i,j];
      e[i,j]:=not e[i,j];
      goto 1;
    end;
    close(dat);
    zapisi_idempotente:=idempotenti;
  end;
function zapisi_idempotente_s_stolpci(s,ime:string):longint;
  label 1,2;
  var
    n,n0,i,j,k:byte;
    e,e2minuse,stevec:array[1..nmax,1..nmax] of boolean;
    vsota:boolean;
    idempotenti:longint;
    dat:text;
  begin
    s:=s+',';
    n:=pos(',',s)-1;
    n0:=length(s) div (n+1);
    for j:=1 to n0 do
    begin
      for i:=1 to n do
      begin
        e[i,j]:=s[1]='1';
        delete(s,1,1);
      end;
      delete(s,1,1);
    end;
    for j:=n0+1 to n do
      for i:=1 to n do
      begin
        stevec[i,j]:=false;
        e[i,j]:=false;
      end;
    for i:=1 to n do
      for j:=1 to n do
      begin
        vsota:=e[i,j];
        for k:=1 to n do
          vsota:=vsota xor (e[i,k] and e[k,j]);
        e2minuse[i,j]:=vsota;
      end;
    assign(dat,ime);
    rewrite(dat);
    idempotenti:=0;
    2:
    for i:=1 to n do
      for j:=1 to n do
        if e2minuse[i,j] then
          goto 1;
    s:='';
    for i:=1 to n do
      for j:=1 to n do
        if e[i,j] then
          s:=s+'1'
        else
          s:=s+'0';
    writeln(dat,s);
    inc(idempotenti);
    1:
    for j:=n0+1 to n do
      for i:=1 to n do
        if stevec[i,j] then
          stevec[i,j]:=false
        else
        begin
          stevec[i,j]:=true;
          for k:=1 to n do
          begin
            if e[k,i] then
              e2minuse[k,j]:=not e2minuse[k,j];
            if e[j,k] then
              e2minuse[i,k]:=not e2minuse[i,k];
          end;
          if i<>j then
            e2minuse[i,j]:=not e2minuse[i,j];
          e[i,j]:=not e[i,j];
          goto 2;
        end;
    close(dat);
    zapisi_idempotente_s_stolpci:=idempotenti;
  end;
function zapisi_idempotente_ver2(n:byte;ime:string):longint;
  label 1,2,3;
  var
    dat:text;
    i,j,k,l,rang,nminusrang:byte;
    ni:array[1..nmax] of byte;
    p:array[1..nmax] of byte;
    a,e0,stevec:array[1..nmax,1..nmax] of boolean;
    s:string;
    idempotenti:longint;
  begin
    assign(dat,ime);
    rewrite(dat);
    idempotenti:=0;
    for rang:=0 to n do
    begin
      nminusrang:=n-rang;
      for i:=1 to rang do
        ni[i]:=i;
      3:
      for i:=1 to rang do
        p[ni[i]]:=i;
      for i:=rang+1 to n do
      begin
        j:=i-rang;
        for k:=1 to rang do
          if ni[k]<=i-rang+k-1 then
            inc(j);
        p[j]:=i;
      end;
      for i:=1 to nminusrang do
        for j:=1 to rang do
          a[i,j]:=false;
      2:
      for i:=1 to rang do
        for j:=1 to nminusrang do
          stevec[i,j]:=false;
      for i:=1 to rang do
      begin
        for j:=1 to rang do
          if i=j then
            e0[i,j]:=true
          else
            e0[i,j]:=false;
        for j:=rang+1 to n do
          e0[i,j]:=false;
      end;
      for i:=rang+1 to n do
      begin
        for j:=1 to rang do
          e0[i,j]:=a[i-rang,j];
        for j:=rang+1 to n do
          e0[i,j]:=false;
      end;
      1:
      s:='';
      for i:=1 to n do
        for j:=1 to n do
          if e0[p[i],p[j]] then
            s:=s+'1'
          else
            s:=s+'0';
      writeln(dat,s);
      inc(idempotenti);
      for i:=1 to rang do
        for j:=1 to nminusrang do
          if stevec[i,j] then
            stevec[i,j]:=false
          else
          begin
            stevec[i,j]:=true;
            e0[i,j+rang]:=not e0[i,j+rang];
            for k:=1 to rang do
              if a[j,k] then
                e0[i,k]:=not e0[i,k];
            for k:=1 to nminusrang do
              for l:=1 to rang do
                if a[k,i] and a[j,l] then
                  e0[k+rang,l]:=not e0[k+rang,l];
            for k:=1 to nminusrang do
              if a[k,i] then
                e0[k+rang,j+rang]:=not e0[k+rang,j+rang];
            goto 1;
          end;
      for j:=1 to rang do
        for i:=1 to ni[j]-j do
          if a[i,j] then
            a[i,j]:=false
          else
          begin
            a[i,j]:=true;
            goto 2;
          end;
      for j:=rang downto 1 do
        if ni[j]<>nminusrang+j then
        begin
          inc(ni[j]);
          for k:=j+1 to rang do
            ni[k]:=ni[j]+k-j;
          goto 3;
        end;
    end;
    close(dat);
    zapisi_idempotente_ver2:=idempotenti;
  end;
procedure hitro_iskanje_ncl3(ime1,ime2:string;pisi:boolean);
  label 1,2,3,4,5;
  var
    dat1,dat2:text;
    n,i,j,k,l,rang,nminusrang:byte;
    ni:array[1..nmax] of byte;
    p:array[1..nmax] of byte;
    a,idem,stevec,x,q:array[1..nmax,1..nmax] of boolean;
    vrstica:array[1..nmax] of boolean;
    vsota:boolean;
    s:string;
    indeks:longint;
  begin
    assign(dat1,ime1);
    assign(dat2,ime2);
    rewrite(dat2);
    close(dat2);
    reset(dat1);
    indeks:=0;
    5:
    while not eof(dat1) do
    begin
      readln(dat1,s);
      inc(indeks);
      n:=round(sqrt(length(s)));
      for i:=1 to n do
        for j:=1 to n do
        begin
          a[i,j]:=s[1]='1';
          delete(s,1,1);
        end;
      rang:=n;
      4:
      nminusrang:=n-rang;
      for i:=1 to rang do
        ni[i]:=i;
      3:
      if pisi then
      begin
        write('Matrika = ',indeks,', Rang = ',rang,', Permutacija = ');
        for i:=1 to rang do
        begin
          write(ni[i]);
          if i<>rang then
            write(',');
        end;
        writeln;
      end;
      for i:=1 to rang do
        p[i]:=ni[i];
      for i:=rang+1 to n do
      begin
        j:=i-rang;
        for k:=1 to rang do
          if ni[k]<=i-rang+k-1 then
            inc(j);
        p[i]:=j;
      end;
      for i:=1 to nminusrang do
        for j:=1 to rang do
          x[i,j]:=false;
      2:
      for i:=1 to rang do
        for j:=1 to nminusrang do
          stevec[i,j]:=false;
      for i:=1 to rang do
      begin
        for j:=1 to rang do
          q[i,j]:=(i=j) xor a[p[i],p[j]];
        for j:=rang+1 to n do
          q[i,j]:=a[p[i],p[j]];
      end;
      for i:=rang+1 to n do
      begin
        for j:=1 to rang do
          q[i,j]:=x[i-rang,j] xor a[p[i],p[j]];
        for j:=rang+1 to n do
          q[i,j]:=a[p[i],p[j]];
      end;
      1:
      for i:=1 to n do
      begin
        for j:=1 to n do
        begin
          vsota:=false;
          for k:=1 to n do
            if q[i,k] and q[k,j] then
              vsota:=not vsota;
          vrstica[j]:=vsota;
        end;
        for j:=1 to n do
        begin
          vsota:=false;
          for k:=1 to n do
            if vrstica[k] and q[k,j] then
              vsota:=not vsota;
          if vsota then
          begin
            for i:=1 to rang do
              for j:=1 to nminusrang do
                if stevec[i,j] then
                  stevec[i,j]:=false
                else
                begin
                  stevec[i,j]:=true;
                  q[i,j+rang]:=not q[i,j+rang];
                  for k:=1 to rang do
                    if x[j,k] then
                      q[i,k]:=not q[i,k];
                  for k:=1 to nminusrang do
                    if x[k,i] then
                      q[k+rang,j+rang]:=not q[k+rang,j+rang];
                  for k:=1 to nminusrang do
                    if x[k,i] then
                      for l:=1 to rang do
                        if x[j,l] then
                          q[k+rang,l]:=not q[k+rang,l];
                  goto 1;
                end;
            for j:=1 to rang do
              for i:=1 to ni[j]-j do
                if x[i,j] then
                  x[i,j]:=false
                else
                begin
                  x[i,j]:=true;
                  goto 2;
                end;
            for j:=rang downto 1 do
              if ni[j]<>nminusrang+j then
              begin
                inc(ni[j]);
                for k:=j+1 to rang do
                  ni[k]:=ni[j]+k-j;
                goto 3;
              end;
            if rang<>0 then
            begin
              dec(rang);
              goto 4;
            end;
            append(dat2);
            writeln(dat2,'');
            close(dat2);
            goto 5;
          end;
        end;
      end;
      for i:=1 to n do
        for j:=1 to n do
          idem[p[i],p[j]]:=q[i,j];
      for i:=1 to n do
        for j:=1 to n do
          idem[i,j]:=idem[i,j] xor a[i,j];
      s:='';
      for i:=1 to n do
        for j:=1 to n do
          if idem[i,j] then
            s:=s+'1'
          else
            s:=s+'0';
      append(dat2);
      writeln(dat2,s);
      close(dat2);
      goto 5;
    end;
    close(dat1);
  end;
function zapisi_idempotente_poddiag(n,r:byte;ime_datoteke:string):longint;
  label 1,2,3;
  var
    i,j,k,j0:byte;
    e,e2minuse,stevec:array[1..nmax,1..nmax] of boolean;
    s:string;
    idempotenti:longint;
    dat:text;
  begin
    for i:=1 to n do
      for j:=1 to n do
      begin
        e[i,j]:=false;
        e2minuse[i,j]:=false;
        stevec[i,j]:=false;
      end;
    assign(dat,ime_datoteke);
    rewrite(dat);
    idempotenti:=0;
    1:
    for i:=1 to n do
      for j:=1 to n do
        if e2minuse[i,j] then
          goto 2;
    s:='';
    for i:=1 to n do
      for j:=1 to n do
        if e[i,j] then
          s:=s+'1'
        else
          s:=s+'0';
    writeln(dat,s);
    inc(idempotenti);
    2:
    for i:=1 to n do
    begin
      if i<=r then
        j0:=1
      else
        j0:=i-r;
      for j:=j0 to n do
        if stevec[i,j] then
          stevec[i,j]:=false
        else
        begin
          stevec[i,j]:=true;
          for k:=1 to n do
          begin
            if e[k,i] then
              e2minuse[k,j]:=not e2minuse[k,j];
            if e[j,k] then
              e2minuse[i,k]:=not e2minuse[i,k];
          end;
          if i<>j then
          e2minuse[i,j]:=not e2minuse[i,j];
          e[i,j]:=not e[i,j];
          goto 1;
        end;
    end;
    close(dat);
    zapisi_idempotente_poddiag:=idempotenti;
  end;
procedure ali_je_ncl(ime1,ime2:string);
  label 1,2,3,4,5;
  var
    dat1,dat2:text;
    s:string;
    n,i,j,k,l,rang,nminusrang:byte;
    ni:array[1..nmax] of byte;
    p,pm1:array[1..nmax] of byte;
    a,papm1,x,e0,stevec,q,q2:array[1..nmax,1..nmax] of boolean;
    vsota:boolean;
  begin
    assign(dat1,ime1);
    assign(dat2,ime2);
    reset(dat1);
    rewrite(dat2);
    while not eof(dat1) do
    begin
      readln(dat1,s);

    n:=round(sqrt(length(s)));
    for i:=1 to n do
      for j:=1 to n do
      begin
        a[i,j]:=s[1]='1';
        delete(s,1,1);
      end;
    rang:=0;
    4:
    nminusrang:=n-rang;
    for i:=1 to rang do
      ni[i]:=i;
    3:
    for j:=1 to rang do
      p[j]:=ni[j];
    for j:=rang+1 to n do
    begin
      i:=j-rang;
      for k:=1 to rang do
        if ni[k]<=j-rang+k-1 then
          inc(i);
      p[j]:=i;
    end;
    for i:=1 to n do
      for j:=1 to n do
        papm1[i,j]:=a[p[i],p[j]];
    for i:=1 to rang do
      for j:=1 to nminusrang do
        x[i,j]:=false;
    2:
    for i:=1 to nminusrang do
      for j:=1 to rang do
        stevec[i,j]:=false;
    for i:=1 to rang do
    begin
      for j:=1 to rang do
        if i=j then
          e0[i,j]:=true
        else
          e0[i,j]:=false;
      for j:=rang+1 to n do
        e0[i,j]:=x[i,j-rang];
    end;
    for i:=rang+1 to n do
      for j:=1 to n do
        e0[i,j]:=false;
    1:
    for i:=1 to n do
      for j:=1 to n do
        q[i,j]:=papm1[i,j] xor e0[i,j];
    for i:=1 to n do
      for j:=1 to n do
      begin
        vsota:=false;
        for k:=1 to n do
          if q[i,k] and q[k,j] then
            vsota:=not vsota;
        q2[i,j]:=vsota;
      end;
    for i:=1 to n do
      for j:=1 to n do
      begin
        vsota:=false;
        for k:=1 to n do
          if q[i,k] and q2[k,j] then
            vsota:=not vsota;
        if vsota then
        begin
          for j:=1 to rang do
            for i:=1 to ni[j]-j do
              if stevec[i,j] then
                stevec[i,j]:=false
              else
              begin
                stevec[i,j]:=true;
                for k:=1 to nminusrang do
                  if e0[k+rang,i+rang] then
                    e0[k+rang,j]:=not e0[k+rang,j];
                for k:=1 to rang do
                  if e0[j,k] then
                    e0[i+rang,k]:=not e0[i+rang,k];
                if x[j,i] then
                  e0[i+rang,j]:=not e0[i+rang,j];
                for k:=1 to nminusrang do
                  if x[j,k] then
                    e0[i+rang,k+rang]:=not e0[i+rang,k+rang];
                for k:=1 to rang do
                  if x[k,i] then
                    e0[k,j]:=not e0[k,j];
                goto 1;
              end;
          for i:=1 to rang do
            for j:=1 to nminusrang do
              if x[i,j] then
                x[i,j]:=false
              else
              begin
                x[i,j]:=true;
                goto 2;
              end;
          for j:=rang downto 1 do
            if ni[j]<>nminusrang+j then
            begin
              inc(ni[j]);
              for k:=j+1 to rang do
                ni[k]:=ni[j]+k-j;
              goto 3;
            end;
          if rang<>n then
          begin
            inc(rang);
            goto 4;
          end;
          s:='';
          goto 5;
        end;
      end;
    for i:=1 to n do
      pm1[p[i]]:=i;
    s:='';
    for i:=1 to n do
      for j:=1 to n do
        if e0[pm1[i],pm1[j]] then
          s:=s+'1'
        else
          s:=s+'0';
    5:
    writeln(dat2,s);
    end;
    close(dat1);
    close(dat2);
  end;
function zapisi_idempotente4(ime1,ime2:string):longint;
  label 1,2,3;
  var
    n,i,j,k,nplus1,vsota:byte;
    e,e2minuse:array[1..nmax+1,1..nmax+1] of byte;
    stevec:array[1..nmax+1,1..nmax+1] of boolean;
    s,z:string;
    idempotenti:longint;
    dat1,dat2:text;
  begin
    assign(dat1,ime1);
    assign(dat2,ime2);
    reset(dat1);
    readln(dat1,s);
    n:=round(sqrt(length(s)));
    nplus1:=n+1;
    close(dat1);
    reset(dat1);
    rewrite(dat2);
    idempotenti:=0;
    while not eof(dat1) do
    begin
      readln(dat1,s);
      for i:=1 to n do
        for j:=1 to n do
        begin
          if s[1]='0' then
            e[i,j]:=0
          else
            e[i,j]:=1;
          delete(s,1,1);
        end;
      for i:=1 to nplus1 do
        e[i,nplus1]:=0;
      for j:=1 to n do
        e[nplus1,j]:=0;
      for i:=1 to n do
        for j:=1 to n do
        begin
          vsota:=0;
          for k:=1 to n do
            vsota:=(vsota+(e[i,k] and e[k,j])) and 3;
          e2minuse[i,j]:=(vsota-e[i,j]) and 3;
        end;
      for i:=1 to nplus1 do
        e2minuse[i,nplus1]:=1;
      for j:=1 to n do
        e2minuse[nplus1,j]:=1;
      for i:=1 to nplus1 do
        for j:=1 to nplus1 do
          stevec[i,j]:=false;
      1:
      i:=1;
      j:=1;
      2:
      if e2minuse[i,j]=0 then
      begin
        if j=n then
        begin
          j:=1;
          inc(i);
        end
        else
          inc(j);
        goto 2;
      end;
      if i=nplus1 then
      begin
        s:='';
        for i:=1 to n do
          for j:=1 to n do
          begin
            str(e[i,j],z);
            s:=s+z;
          end;
        writeln(dat2,s);
        inc(idempotenti);
      end;
      i:=1;
      j:=1;
      3:
      if stevec[i,j] then
      begin
        stevec[i,j]:=false;
        if j=n then
        begin
          j:=1;
          inc(i);
        end
        else
          inc(j);
        goto 3;
      end;
      if i<>nplus1 then
      begin
        stevec[i,j]:=true;
        for k:=1 to n do
        begin
          e2minuse[k,j]:=(e2minuse[k,j]+(e[k,i] shl 1)) and 3;
          e2minuse[i,k]:=(e2minuse[i,k]+(e[j,k] shl 1)) and 3;
        end;
        e2minuse[i,j]:=(e2minuse[i,j]+2) and 3;
        e[i,j]:=e[i,j] xor 2;
        goto 1;
      end;
    end;
    close(dat1);
    close(dat2);
    zapisi_idempotente4:=idempotenti;
  end;
function zapisi_frobeniusove_matrike(velikosti,ime:string):longint;
  label 1;
  var
    k,n,r:byte;
    vel:array[1..nmax] of byte;
    i,j,koda:integer;
    stevec:array[1..nmax+1] of boolean;
    matrika:array[1..nmax,1..nmax] of boolean;
    stevilo:longint;
    dat:text;
    s:string;
  begin
    k:=0;
    n:=0;
    velikosti:=velikosti+',';
    while velikosti<>'' do
    begin
      inc(k);
      i:=pos(',',velikosti);
      val(copy(velikosti,1,i-1),vel[k],koda);
      delete(velikosti,1,i);
      n:=n+vel[k];
    end;
    for i:=1 to n do
      for j:=1 to n do
        matrika[i,j]:=false;
    r:=0;
    for i:=1 to k do
    begin
      for j:=1 to vel[i]-1 do
        matrika[r+j+1,r+j]:=true;
      r:=r+vel[i];
    end;
    for i:=1 to n+1 do
      stevec[i]:=false;
    assign(dat,ime);
    rewrite(dat);
    stevilo:=0;
    1:
    r:=0;
    for i:=1 to k do
    begin
      for j:=1 to vel[i] do
        matrika[r+j,r+vel[i]]:=stevec[r+j];
      r:=r+vel[i];
    end;
    s:='';
    for i:=1 to n do
      for j:=1 to n do
        if matrika[i,j] then
          s:=s+'1'
        else
          s:=s+'0';
    writeln(dat,s);
    inc(stevilo);
    i:=1;
    while stevec[i] do
    begin
      stevec[i]:=false;
      inc(i);
    end;
    stevec[i]:=true;
    if i<>n+1 then
      goto 1;
    close(dat);
    zapisi_frobeniusove_matrike:=stevilo;
  end;
function zapisi_mod_frobeniusove_matrike(velikosti,ime:string):longint;
  label 1;
  var
    k,n,r:byte;
    vel:array[1..nmax] of byte;
    i,j,koda:integer;
    stevec:array[1..nmax+1] of boolean;
    matrika:array[1..nmax,1..nmax] of boolean;
    stevilo:longint;
    dat:text;
    s:string;
  begin
    k:=0;
    n:=0;
    velikosti:=velikosti+',';
    while velikosti<>'' do
    begin
      inc(k);
      i:=pos(',',velikosti);
      val(copy(velikosti,1,i-1),vel[k],koda);
      delete(velikosti,1,i);
      n:=n+vel[k];
    end;
    for i:=1 to n do
      for j:=1 to n do
        matrika[i,j]:=false;
    r:=0;
    for i:=1 to k do
    begin
      for j:=1 to vel[i]-1 do
      begin
        matrika[r+j+1,r+j]:=true;
        if j mod 2=1 then
          matrika[r+j,r+j]:=true;
      end;
      r:=r+vel[i];
    end;
    for i:=1 to n+1 do
      stevec[i]:=false;
    assign(dat,ime);
    rewrite(dat);
    stevilo:=0;
    1:
    r:=0;
    for i:=1 to k do
    begin
      for j:=1 to vel[i] do
        matrika[r+j,r+vel[i]]:=stevec[r+j];
      r:=r+vel[i];
    end;
    s:='';
    for i:=1 to n do
      for j:=1 to n do
        if matrika[i,j] then
          s:=s+'1'
        else
          s:=s+'0';
    writeln(dat,s);
    inc(stevilo);
    i:=1;
    while stevec[i] do
    begin
      stevec[i]:=false;
      inc(i);
    end;
    stevec[i]:=true;
    if i<>n+1 then
      goto 1;
    close(dat);
    zapisi_mod_frobeniusove_matrike:=stevilo;
  end;
function zapisi_problematicne_mod_frobeniusove_matrike(velikosti,ime:string):longint;
  label 1;
  const
    problematicna1='1001100001110011';
    problematicna2='1001100001100011';
  var
    k,n,r:byte;
    vel:array[1..nmax] of byte;
    i,j,koda:integer;
    stevec:array[1..nmax] of boolean;
    matrika:array[1..nmax,1..nmax] of boolean;
    stevilo:longint;
    dat:text;
    s:string;
  begin
    k:=0;
    n:=0;
    velikosti:=velikosti+',';
    while velikosti<>'' do
    begin
      inc(k);
      i:=pos(',',velikosti);
      val(copy(velikosti,1,i-1),vel[k],koda);
      delete(velikosti,1,i);
      n:=n+vel[k];
    end;
    for i:=1 to n+4 do
      for j:=1 to n+4 do
        matrika[i,j]:=false;
    r:=0;
    for i:=1 to k do
    begin
      for j:=1 to vel[i]-1 do
      begin
        matrika[r+j+1,r+j]:=true;
        if j mod 2=1 then
          matrika[r+j,r+j]:=true;
      end;
      r:=r+vel[i];
    end;
    for i:=1 to n do
      stevec[i]:=false;
    assign(dat,ime);
    rewrite(dat);
    stevilo:=0;
    1:
    r:=0;
    for i:=1 to k do
    begin
      for j:=1 to vel[i] do
        matrika[r+j,r+vel[i]]:=stevec[r+j];
      r:=r+vel[i];
    end;
    s:='';
    for i:=1 to n do
    begin
      for j:=1 to n do
        if matrika[i,j] then
          s:=s+'1'
        else
          s:=s+'0';
      s:=s+'0000';
    end;
    for i:=1 to 4 do
    begin
      for j:=1 to n do
        s:=s+'0';
      s:=s+copy(problematicna1,4*i-3,4);
    end;
    writeln(dat,s);
    inc(stevilo);
    s:='';
    for i:=1 to n do
    begin
      for j:=1 to n do
        if matrika[i,j] then
          s:=s+'1'
        else
          s:=s+'0';
      s:=s+'0000';
    end;
    for i:=1 to 4 do
    begin
      for j:=1 to n do
        s:=s+'0';
      s:=s+copy(problematicna2,4*i-3,4);
    end;
    writeln(dat,s);
    inc(stevilo);
    for i:=1 to n do
      if stevec[i] then
        stevec[i]:=false
      else
      begin
        stevec[i]:=true;
        goto 1;
      end;
    close(dat);
    zapisi_problematicne_mod_frobeniusove_matrike:=stevilo;
  end;
procedure poisci_stevilo_ncl_dekompozicij(ime1,ime2,ime3:string);
  label 1;
  var
    matrika1,matrika2,matrika3:array[1..nmax,1..nmax] of boolean;
    dat1,dat2,dat3:text;
    n,i,j,k:byte;
    s:string;
    stevilo:longint;
    vsota:boolean;
  begin
    assign(dat1,ime1);
    assign(dat2,ime2);
    assign(dat3,ime3);
    reset(dat1);
    rewrite(dat3);
    while not eof(dat1) do
    begin
      readln(dat1,s);
      n:=round(sqrt(length(s)));
      for i:=1 to n do
        for j:=1 to n do
        begin
          if s[1]='1' then
            matrika1[i,j]:=true
          else
            matrika1[i,j]:=false;
          delete(s,1,1);
        end;
      reset(dat2);
      stevilo:=0;
      while not eof(dat2) do
      begin
        readln(dat2,s);
        for i:=1 to n do
          for j:=1 to n do
          begin
            if s[1]='1' then
              matrika2[i,j]:=not matrika1[i,j]
            else
              matrika2[i,j]:=matrika1[i,j];
            delete(s,1,1);
          end;
        for i:=1 to n do
          for j:=1 to n do
          begin
            vsota:=false;
            for k:=1 to n do
              vsota:=vsota xor (matrika2[i,k] and matrika2[k,j]);
            matrika3[i,j]:=vsota;
          end;
        i:=1;
        j:=1;
        1:
        vsota:=false;
        for k:=1 to n do
          vsota:=vsota xor (matrika2[i,k] and matrika3[k,j]);
        if not vsota then
        begin
          if j=n then
          begin
            j:=1;
            inc(i);
          end
          else
            inc(j);
          if i<>n+1 then
            goto 1;
        end;
        if i=n+1 then
          inc(stevilo);
      end;
      close(dat2);
      writeln(dat3,stevilo);
    end;
    close(dat1);
    close(dat3);
  end;
function preveri_ncl3(ime1,ime2:string):boolean;
  label 1,2;
  var
    matrika1,matrika2:array[1..nmax,1..nmax] of boolean;
    vrstica:array[1..nmax] of boolean;
    dat1,dat2:text;
    n,i,j,k:byte;
    s:string;
    vsota:boolean;
    u:boolean;
  begin
    assign(dat1,ime1);
    assign(dat2,ime2);
    reset(dat1);
    while not eof(dat1) do
    begin
      readln(dat1,s);
      n:=round(sqrt(length(s)));
      for i:=1 to n do
        for j:=1 to n do
        begin
          if s[1]='1' then
            matrika1[i,j]:=true
          else
            matrika1[i,j]:=false;
          delete(s,1,1);
        end;
      reset(dat2);
      u:=false;
      while not u and not eof(dat2) do
      begin
        readln(dat2,s);
        for i:=1 to n do
          for j:=1 to n do
          begin
            if s[1]='1' then
              matrika2[i,j]:=not matrika1[i,j]
            else
              matrika2[i,j]:=matrika1[i,j];
            delete(s,1,1);
          end;
        for i:=1 to n do
        begin
          for j:=1 to n do
          begin
            vsota:=false;
            for k:=1 to n do
              if matrika2[i,k] and matrika2[k,j] then
                vsota:=not vsota;
            vrstica[j]:=vsota;
          end;
          for j:=1 to n do
          begin
            vsota:=false;
            for k:=1 to n do
              if vrstica[k] and matrika2[k,j] then
                vsota:=not vsota;
            if vsota then
              goto 1;
          end;
        end;
        u:=true;
        1:
      end;
      close(dat2);
      if not u then
      begin
        close(dat1);
        preveri_ncl3:=false;
        goto 2;
      end;
    end;
    close(dat1);
    preveri_ncl3:=true;
    2:
  end;
function najdi_idempotente(s,ime1,ime2:string):longint;
  label 1;
  var
    n,i,j,k:byte;
    z:string;
    a,b,c:array[1..nmax,1..nmax] of boolean;
    vsota:boolean;
    dat1,dat2:text;
    idempotenti:longint;
  begin
    n:=round(sqrt(length(s)));
    for i:=1 to n do
      for j:=1 to n do
      begin
        if s[1]='0' then
          a[i,j]:=false
        else
          a[i,j]:=true;
        delete(s,1,1);
      end;
    idempotenti:=0;
    assign(dat1,ime1);
    assign(dat2,ime2);
    reset(dat1);
    rewrite(dat2);
    while not eof(dat1) do
    begin
      readln(dat1,s);
      z:=s;
      for i:=1 to n do
        for j:=1 to n do
        begin
          if s[1]='0' then
            b[i,j]:=a[i,j]
          else
            b[i,j]:=not a[i,j];
          delete(s,1,1);
        end;
      for i:=1 to n do
        for j:=1 to n do
        begin
          vsota:=false;
          for k:=1 to n do
            vsota:=vsota xor (b[i,k] and b[k,j]);
          c[i,j]:=vsota;
        end;
      for i:=1 to n do
        for j:=1 to n do
        begin
          vsota:=false;
          for k:=1 to n do
            vsota:=vsota xor (c[i,k] and c[k,j]);
          if vsota then
            goto 1;
        end;
      writeln(dat2,z);
      inc(idempotenti);
      1:
    end;
    close(dat1);
    close(dat2);
    najdi_idempotente:=idempotenti;
  end;
procedure poisci_stevilo_ncl_dekompozicij4(s,ime1,ime2:string;r:byte);
  label 1,2;
  var
    n,i,j,k,t,vsota:byte;
    a:array[1..nmax,1..nmax] of byte;
    q:array[1..nmax,1..nmax,1..maxpower] of byte;
    stevec:array[1..nmax,1..nmax] of boolean;
    dat1,dat2:text;
    z:string;
    idempotenti:longint;
  begin
    n:=round(sqrt(length(s)));
    for i:=1 to n do
      for j:=1 to n do
      begin
        if s[1]='0' then
          a[i,j]:=0
        else
          a[i,j]:=1;
        delete(s,1,1);
        stevec[i,j]:=false;
      end;
    assign(dat1,ime1);
    assign(dat2,ime2);
    rewrite(dat2);
    2:
    s:='';
    for i:=1 to n do
      for j:=1 to n do
      begin
        str(a[i,j] and 3,z);
        s:=s+z;
      end;
    write(dat2,s+': ');
    idempotenti:=0;
    reset(dat1);
    while not eof(dat1) do
    begin
      readln(dat1,s);
      for i:=1 to n do
        for j:=1 to n do
        begin
          z:=s[1];
          if z='0' then
            q[i,j,1]:=a[i,j]
          else if z='1' then
            q[i,j,1]:=a[i,j]-1
          else if z='2' then
            q[i,j,1]:=a[i,j]-2
          else
            q[i,j,1]:=a[i,j]-3;
          delete(s,1,1);
        end;
      for t:=2 to r do
        for i:=1 to n do
          for j:=1 to n do
          begin
            vsota:=0;
            for k:=1 to n do
              vsota:=vsota+q[i,k,t-1]*q[k,j,1];
            q[i,j,t]:=vsota;
          end;
      for i:=1 to n do
        for j:=1 to n do
          if (q[i,j,r] and 3)<>0 then
            goto 1;
      inc(idempotenti);
      1:
    end;
    close(dat1);
    writeln(dat2,idempotenti);
    for i:=1 to n do
      for j:=1 to n do
        if stevec[i,j] then
          stevec[i,j]:=false
        else
        begin
          stevec[i,j]:=true;
          a[i,j]:=a[i,j]+2;
          goto 2;
        end;
    close(dat2);
  end;
function najdi_idempotente4(s,ime1,ime2:string;r:byte):longint;
  label 1;
  var
    n,i,j,k,t,vsota:byte;
    a:array[1..nmax,1..nmax] of byte;
    qn:array[1..nmax,1..nmax,1..8] of byte;
    dat1,dat2:text;
    s1,z:string;
    idempotenti:longint;
  begin
    n:=round(sqrt(length(s)));
    for i:=1 to n do
      for j:=1 to n do
      begin
        z:=s[1];
        if z='0' then
          a[i,j]:=0
        else if z='1' then
          a[i,j]:=1
        else if z='2' then
          a[i,j]:=2
        else
          a[i,j]:=3;
        delete(s,1,1);
      end;
    assign(dat1,ime1);
    assign(dat2,ime2);
    reset(dat1);
    rewrite(dat2);
    idempotenti:=0;
    while not eof(dat1) do
    begin
      readln(dat1,s1);
      s:=s1;
      for i:=1 to n do
        for j:=1 to n do
        begin
          z:=s[1];
          if z='0' then
            qn[i,j,1]:=a[i,j]
          else if z='1' then
            qn[i,j,1]:=a[i,j]-1
          else if z='2' then
            qn[i,j,1]:=a[i,j]-2
          else
            qn[i,j,1]:=a[i,j]-3;
          delete(s,1,1);
        end;
      for t:=2 to r do
        for i:=1 to n do
          for j:=1 to n do
          begin
            vsota:=0;
            for k:=1 to n do
              vsota:=vsota+qn[i,k,t-1]*qn[k,j,1];
            qn[i,j,t]:=vsota;
          end;
      for i:=1 to n do
        for j:=1 to n do
          if (qn[i,j,r] and 3)<>0 then
            goto 1;
      writeln(dat2,s1);
      inc(idempotenti);
      1:
    end;
    close(dat1);
    close(dat2);
    najdi_idempotente4:=idempotenti;
  end;
function primerjaj_datoteki(ime1,ime2:string):boolean;
  label 1,2;
  var
    s,z:string;
    dat1,dat2:text;
  begin
    assign(dat1,ime1);
    assign(dat2,ime2);
    reset(dat1);
    while not eof(dat1) do
    begin
      readln(dat1,s);
      reset(dat2);
      while not eof(dat2) do
      begin
        readln(dat2,z);
        if s=z then
        begin
          close(dat2);
          goto 1;
        end;
      end;
      close(dat2);
      close(dat1);
      primerjaj_datoteki:=false;
      goto 2;
      1:
    end;
    close(dat1);
    primerjaj_datoteki:=true;
    2:
  end;
procedure predstavitev(ime1,ime2,ime3:string);
  const
    znak='1';
  var
    dat1,dat2,dat3:text;
    s,z:string;
    n,i,j:byte;
    a,e,q:array[1..nmax,1..nmax] of boolean;
  begin
    assign(dat1,ime1);
    assign(dat2,ime2);
    assign(dat3,ime3);
    reset(dat1);
    reset(dat2);
    rewrite(dat3);
    while not eof(dat1) do
    begin
      readln(dat1,s);
      readln(dat2,z);
      n:=round(sqrt(length(s)));
      for i:=1 to n do
        for j:=1 to n do
        begin
          a[i,j]:=s[1]='1';
          e[i,j]:=z[1]='1';
          q[i,j]:=a[i,j] xor e[i,j];
          delete(s,1,1);
          delete(z,1,1);
        end;
      for i:=1 to n do
      begin
        s:='[';
        for j:=1 to n-1 do
          if a[i,j] then
            s:=s+znak+' '
          else
            s:=s+'  ';
        if a[i,n] then
          s:=s+znak
        else
          s:=s+' ';
        if i=(n+1) div 2 then
          s:=s+']   =   ['
        else
          s:=s+']       [';
        for j:=1 to n-1 do
          if e[i,j] then
            s:=s+znak+' '
          else
            s:=s+'  ';
        if e[i,n] then
          s:=s+znak
        else
          s:=s+' ';
        if i=(n+1) div 2 then
          s:=s+']   +   ['
        else
          s:=s+']       [';
        for j:=1 to n-1 do
          if q[i,j] then
            s:=s+znak+' '
          else
            s:=s+'  ';
        if q[i,n] then
          s:=s+znak+']'
        else
          s:=s+' ]';
        writeln(dat3,s);
      end;
      writeln(dat3,'');
    end;
    close(dat1);
    close(dat2);
    close(dat3);
  end;

  label 1;
  var
    n,r,ver:byte;
    s,ime,ime1,ime2,ime3,velikosti:string;
    stevilo:longint;
    u:boolean;
  begin
    1:
    clrscr;
    writeln(' 1 = Operacije modulo 2');
    writeln(' 2 = Operacije modulo 4');
    writeln(' 3 = Splosne operacije');
    writeln(' 4 = Izhod');
    readln(n);
    clrscr;
    if n=1 then
    begin
      writeln(' 1 = Zapisi vse idempotente');
      writeln(' 2 = Zapisi vse Frobeniusove matrike');
      writeln(' 3 = Preveri obstoj nil-clean dekompozicij indeksa 3');
      writeln(' 4 = Poisci stevilo nil-clean dekompozicij indeksa 3');
      writeln(' 5 = Poisci vse idempotente v nil-clean dekompozicijah indeksa 4');
      writeln(' 6 = Poisci idempotent v kaki dekompoziciji indeksa 3 za dane matrike');
      writeln(' 7 = Pripravi uporabniku prijazno predstavitev dekompozicij');
      writeln(' 8 = Zapisi vse idempotente, ki so nic pod r-to poddiagonalo');
      writeln(' 9 = Zapisi vse blocno zgornje trikotne idempotente z danimi bloki');
      writeln('10 = Zapisi vse modificirane Frobeniusove matrike');
      writeln('11 = Za dano matriko poisci vse NCL idempotente iz datoteke');
      writeln('12 = Zapisi modificirane Frobeniusove matrike s problematicnimi 4x4 bloki');
      writeln('13 = Zapisi vse idempotente z danimi stolpci');
      writeln('14 = Zapisi vse blocno zgornje trikotne idempotente z danima blokoma');
      writeln('15 = Hitro iskanje idempotenta v NCL dekompoziciji reda 3 za dane matrike');
      writeln('16 = Glavni meni');
      readln(n);
      clrscr;
      if n=1 then
      begin
        writeln('Program bo poiskal vse idempotentne matrike dane velikosti.');
        write('Vnesi velikost matrik: ');
        readln(n);
        write('Vnesi ime ciljne datoteke: ');
        readln(ime);
        write('Vnesi verzijo algoritma (1 = preveri vse, 2 = fast): ');
        readln(ver);
        write('Racunam ...');
        if ver=1 then
          stevilo:=zapisi_idempotente(n,ime)
        else
          stevilo:=zapisi_idempotente_ver2(n,ime);
        writeln(' Koncano.');
        writeln('Stevilo najdenih idempotentov: ',stevilo);
      end
      else if n=2 then
      begin
        writeln('Program bo poiskal vse Frobeniusove matrike z bloki danih velikosti.');
        write('Vnesi velikosti blokov (velikosti naj bodo locene z vejicami): ');
        readln(velikosti);
        write('Vnesi ime ciljne datoteke: ');
        readln(ime);
        write('Racunam ...');
        stevilo:=zapisi_frobeniusove_matrike(velikosti,ime);
        writeln(' Koncano.');
        writeln('Stevilo najdenih matrik: ',stevilo);
      end
      else if n=3 then
      begin
        writeln('Program bo preveril, ali za vsako matriko X iz datoteke A obstaja matrika Y iz');
        writeln('datoteke B, da velja (X-Y)^3=0.');
        write('Vnesi ime datoteke A: ');
        readln(ime1);
        write('Vnesi ime datoteke B: ');
        readln(ime2);
        write('Racunam ...');
        u:=preveri_ncl3(ime1,ime2);
        writeln(' Koncano.');
        if u then
          writeln('Odgovor je DA.')
        else
          writeln('Odgovor je NE.');
      end
      else if n=4 then
      begin
        writeln('Program bo za vsako matriko X iz datoteke A poiskal stevilo matrik Y iz');
        writeln('datoteke B, za katere velja (X-Y)^3=0. Dobljena stevila bo zapisal v datoteko');
        writeln('C.');
        write('Vnesi ime datoteke A: ');
        readln(ime1);
        write('Vnesi ime datoteke B: ');
        readln(ime2);
        write('Vnesi ime datoteke C: ');
        readln(ime3);
        write('Racunam ...');
        poisci_stevilo_ncl_dekompozicij(ime1,ime2,ime3);
        writeln(' Koncano.');
      end
      else if n=5 then
      begin
        writeln('Program bo za dano matriko X poiskal vse matrike Y iz datoteke A, za katere');
        writeln('velja (X-Y)^4=0. Matrike bo zapisal v datoteko B.');
        write('Vnesi matriko X: ');
        readln(s);
        write('Vnesi ime datoteke A: ');
        readln(ime1);
        write('Vnesi ime datoteke B: ');
        readln(ime2);
        write('Racunam ...');
        stevilo:=najdi_idempotente(s,ime1,ime2);
        writeln(' Koncano.');
        writeln('Stevilo najdenih idempotentov: ',stevilo);
      end
      else if n=6 then
      begin
        writeln('Program bo za vsako matriko X iz datoteke A poiskal idempotent E z lastnostjo');
        writeln('(X-E)^3=0 (ce obstaja). Idempotente bo zapisal v datoteko B.');
        write('Vnesi ime datoteke A: ');
        readln(ime1);
        write('Vnesi ime datoteke B: ');
        readln(ime2);
        write('Racunam ...');
        ali_je_ncl(ime1,ime2);
        writeln(' Koncano.');
      end
      else if n=7 then
      begin
        writeln('Program bo po vrsti bral matrike iz datotek A in B; ce sta X in Y ti matriki,');
        writeln('bo program izracunal Z=X-Y in v datoteko C zapisal uporabniku prijazno');
        writeln('prestavitev enacbe X=Y+Z.');
        write('Vnesi ime datoteke A: ');
        readln(ime1);
        write('Vnesi ime datoteke B: ');
        readln(ime2);
        write('Vnesi ime datoteke C: ');
        readln(ime3);
        write('Zapisujem ...');
        predstavitev(ime1,ime2,ime3);
        writeln(' Koncano.');
      end
      else if n=8 then
      begin
        writeln('Program bo v datoteko A zapisal vse idempotente, ki so 0 pod r-to poddiagonalo.');
        write('Vnesi velikost matrik: ');
        readln(n);
        write('Vnesi r: ');
        readln(r);
        write('Vnesi ime datoteke A: ');
        readln(ime);
        write('Racunam ...');
        stevilo:=zapisi_idempotente_poddiag(n,r,ime);
        writeln(' Koncano.');
        writeln('Stevilo najdenih idempotentov: ',stevilo);
      end
      else if n=9 then
      begin
        writeln('Program bo poiskal in zapisal v datoteko A vse blocno zgornje trikotne');
        writeln('idempotente z vnaprej podanimi diagonalnimi bloki.');
        write('Vnesi diagonalne bloke (locene z vejicami): ');
        readln(s);
        write('Vnesi ime datoteke A: ');
        readln(ime);
        write('Racunam ...');
        stevilo:=zapisi_bzt_idempotente(s,ime);
        writeln(' Koncano.');
        writeln('Stevilo najdenih idempotentov: ',stevilo);
      end
      else if n=10 then
      begin
        writeln('Program bo v datoteko A zapisal vse modificirane Frobeniusove matrike z bloki');
        writeln('danih velikosti.');
        write('Vnesi velikosti blokov (locene z vejicami): ');
        readln(s);
        write('Vnesi ime datoteke A: ');
        readln(ime);
        write('Racunam ...');
        stevilo:=zapisi_mod_frobeniusove_matrike(s,ime);
        writeln(' Koncano.');
        writeln('Stevilo najdenih matrik: ',stevilo);
      end
      else if n=11 then
      begin
        writeln('Program bo za dano matriko X poiskal vse matrike Y iz datoteke A, ki zadoscajo');
        writeln('(X-Y)^3=0. Matrike Y bo zapisal v datoteko B.');
        write('Vnesi matriko X: ');
        readln(s);
        write('Vnesi ime datoteke A: ');
        readln(ime1);
        write('Vnesi ime datoteke B: ');
        readln(ime2);
        write('Racunam ...');
        stevilo:=izlusci_idempotente(s,ime1,ime2);
        writeln(' Koncano.');
        writeln('Stevilo najdenih matrik: ',stevilo);
      end
      else if n=12 then
      begin
        writeln('Program bo v datoteko A zapisal vse modificirane Frobeniusove matrike, ki');
        writeln('vsebujejo kakega od problematicnih 4x4 blokov. Velikosti preostalih blokov');
        writeln('poda uporabnik.');
        write('Vnesi velikosti ostalih blokov: ');
        readln(velikosti);
        write('Vnesi ime datoteke A: ');
        readln(ime);
        write('Racunam ...');
        stevilo:=zapisi_problematicne_mod_frobeniusove_matrike(velikosti,ime);
        writeln(' Koncano.');
        writeln('Stevilo najdenih matrik: ',stevilo);
      end
      else if n=13 then
      begin
        writeln('Program bo v datoteko zapisal vse idempotente z danimi prvimi nekaj stolpci.');
        write('Vnesi prve stolpce (locene z vejicami): ');
        readln(s);
        write('Vnesi ime datoteke: ');
        readln(ime);
        write('Racunam ...');
        stevilo:=zapisi_idempotente_s_stolpci(s,ime);
        writeln(' Koncano.');
        writeln('Stevilo najdenih idempotentov: ',stevilo);
      end
      else if n=14 then
      begin
        writeln('Program bo za vsak idempotent X iz datoteke A in idempotent Y iz datoteke B');
        writeln('poiskal vse blocno zgornje trikotne idempotente z diagonalnima blokoma X in Y.');
        writeln('Idempotente bo zapisal v datoteko C.');
        write('Vnesi ime datoteke A: ');
        readln(ime1);
        write('Vnesi ime datoteke B: ');
        readln(ime2);
        write('Vnesi ime datoteke C: ');
        readln(ime3);
        write('Racunam ...');
        stevilo:=blocno_zgornje_trikotni_idempotenti(ime1,ime2,ime3);
        writeln(' Koncano.');
        writeln('Stevilo najdenih idempotentov: ',stevilo);
      end
      else if n=15 then
      begin
        writeln('Program bo za vsako matriko X iz datoteke A poiskal kak idempotent E, za');
        writeln('katerega velja (X-E)^3=0. Idempotente bo zapisal v datoteko B.');
        write('Vnesi ime datoteke A: ');
        readln(ime1);
        write('Vnesi ime datoteke B: ');
        readln(ime2);
        write('Racunam ...');
        hitro_iskanje_ncl3(ime1,ime2,true);
        writeln(' Koncano.');
      end
      else
        goto 1;
    end
    else if n=2 then
    begin
      writeln(' 1 = Poisci vse idempotente modulo 4 iz danih idempotentov modulo 2');
      writeln(' 2 = Poisci stevilo idempotentov v NCL dek. indeksa n dane matrike modulo 2');
      writeln(' 3 = Poisci vse idempotente v NCL dek. indeksa n dane matrike');
      writeln(' 4 = Glavni meni');
      readln(n);
      clrscr;
      if n=1 then
      begin
        writeln('Program bo poiskal vse idempotente modulo 4, ki so modulo 2 enaki kakemu izmed');
        writeln('idempotentov, ki se nahajajo v datoteki A. Idempotente bo zapisal v datoteko B.');
        write('Vnesi ime datoteke A: ');
        readln(ime1);
        write('Vnesi ime datoteke B: ');
        readln(ime2);
        write('Racunam ...');
        stevilo:=zapisi_idempotente4(ime1,ime2);
        writeln(' Koncano.');
        writeln('Stevilo najdenih idempotentov: ',stevilo);
      end
      else if n=2 then
      begin
        writeln('Program bo za dano matriko X modulo 2 poiskal vse njene predstavnike modulo 4.');
        writeln('Nato bo za vsakega predstavnika Y poiskal stevilo matrik Z iz datoteke A,');
        writeln('za katere velja (Y-Z)^n=0. Dobljena stevila bo zapisal v datoteko B.');
        write('Vnesi matriko X: ');
        readln(s);
        write('Vnesi ime datoteke A: ');
        readln(ime1);
        write('Vnesi ime datoteke B: ');
        readln(ime2);
        write('Vnesi n: ');
        readln(n);
        write('Racunam ...');
        poisci_stevilo_ncl_dekompozicij4(s,ime1,ime2,n);
        writeln(' Koncano.');
      end
      else if n=3 then
      begin
        writeln('Program bo za dano matriko X poiskal vse matrike Y iz datoteke A, za katere');
        writeln('velja (X-Y)^n=0. Matrike Y bo zapisal v datoteko B.');
        write('Vnesi matriko X: ');
        readln(s);
        write('Vnesi ime datoteke A: ');
        readln(ime1);
        write('Vnesi ime datoteke B: ');
        readln(ime2);
        write('Vnesi n: ');
        readln(n);
        write('Racunam ...');
        stevilo:=najdi_idempotente4(s,ime1,ime2,n);
        writeln(' Koncano.');
        writeln('Stevilo najdenih idempotentov: ',stevilo);
      end
      else
        goto 1;
    end
    else if n=3 then
    begin
      writeln(' 1 = Primerjaj datoteki');
      writeln(' 2 = Glavni meni');
      readln(n);
      clrscr;
      if n=1 then
      begin
        writeln('Program bo preveril, ali se vsaka vrstica datoteke A nahaja tudi v datoteki B.');
        write('Vnesi ime datoteke A: ');
        readln(ime1);
        write('Vnesi ime datoteke B: ');
        readln(ime2);
        write('Primerjam ...');
        u:=primerjaj_datoteki(ime1,ime2);
        writeln(' Koncano.');
        if u then
          writeln('Odgovor je DA.')
        else
          writeln('Odgovor je NE.');
      end
      else
        goto 1;
    end
    else if n=4 then
      exit
    else
      goto 1;
    writeln('Pritisni Enter za vrnitev na glavni meni.');
    readln;
    goto 1;
  end.