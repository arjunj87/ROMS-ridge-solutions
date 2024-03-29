      implicit none
      integer max_length, max_filename
      parameter(max_length=72, max_filename=32)
      character file_in*(max_filename), file_out*(max_filename)
      integer iin, iout, last_arg, iargc, itype
 
      character(len=32) varname(1024), svarname(1024),  vardim(1024)
      integer ivar, nvars, namelen(1024), dimlen(1024), varstat(1024)
      data nvars/0/ namelen/1024*0/ dimlen/1024*0/ varstat/1024*0/
 
      character(len=max_length) string, scratch
      character quote
      logical not_end_of_file, lswtch, directive
      integer line, length, ierr, i,j, istrt, ir, ig, nwords,
     &      is(32),ie(32), io,ic
     &   ,case_fold, ibrack
 
      quote=char(39)
      case_fold=ichar('a')-ichar('A')
      ibrack=0
 
      file_in='                                '
      file_out='                                '
      last_arg=iargc()
      if (last_arg.gt.0) then
        iin=11
        call getarg(1,file_in)
        open(iin, file=file_in, form='formatted', status='old')
      else
        iin=5       !<-- read from standard input
      endif
      if (last_arg.eq.2) then
        iout=12
        call getarg(2,file_out)
        open(iout, file=file_out, form='formatted', status='unknown')
      else
        iout=6        !--> write to standard output
      endif
 
      not_end_of_file=.true.
      line=0
 
      do while (not_end_of_file)
        line=line+1
        length=0
        string=' '
        read(iin,'(A)',iostat=ierr,end=1) string
        goto 2
   1    not_end_of_file=.false.
   2    length=max_length
        do while (length.gt.0 .and. string(length:length).eq.' ')
          length=length-1
        enddo
        directive=.false.
        if (string(1:1).eq.'!' .or. string(1:1).eq.'C' .or.
     &                              string(1:1).eq.'c') then
          do i=2,5
            if (string(i:i).eq.'$') directive=.true.
          enddo
          if (.not.directive) length=0
        endif
        istrt=1
        do while (string(istrt:istrt).eq.' ' .and. istrt.lt.length)
          istrt=istrt+1
        enddo
        if (.not.directive .and. string(istrt:istrt).eq.'!') length=0
        if (istrt.le.length .and. string(istrt:istrt).ne.'#') then
          i=0
          lswtch=.false.
          do while (i.lt.length)
            i=i+1
            if (string(i:i).eq.quote) then
              lswtch=.not.lswtch
            elseif (string(i:i).eq.'!' .and. .not.lswtch) then
              length=i-1
            endif
          enddo
!
! Extract names of declared variables: to do so create a temporal
! string "scratch", which is equal to the original, except that
! upper/lower cases, commas, expressions in brackets (possibly array
! dimensions) are all suppressed. Once this is done, classify line
! to be either FORTRAN type specifier, common block, or parameter
! declaration. Move starting index of the string to prevent it from
! treating FORTRAN reserved names as variable names. Note that type
! of the line "itype" is retained from the previous line, if the
! present one is a continuation line.
!
          scratch(1:istrt)=' '
          ibrack=0
          do i=istrt,length
            if (string(i:i).eq.'(') then
              scratch(i:i)=' '
              ibrack=ibrack+1
            elseif (string(i:i).eq.')') then
              scratch(i:i)=' '
              ibrack=ibrack-1
            elseif (ibrack.gt.0) then
              scratch(i:i)=' '
            elseif (string(i:i).eq.',') then
              scratch(i:i)=' '
            elseif (string(i:i).ge.'A' .and. string(i:i).le.'Z') then
              scratch(i:i)=char(ichar(string(i:i))+case_fold)
            else
              scratch(i:i)=string(i:i)
            endif
          enddo
          scratch(length+1:max_length)=' '
 
          if (string(6:6).eq.' ') then
            if (scratch(istrt:istrt+3).eq.'real') then
              itype=1
              istrt=istrt+4
            elseif (scratch(istrt:istrt+6).eq.'integer') then
              itype=2
              istrt=istrt+7
            elseif (scratch(istrt:istrt+8).eq.'character') then
              itype=3
              istrt=istrt+9
            elseif (scratch(istrt:istrt+8).eq.'parameter') then
              itype=-1
              istrt=istrt+9
            elseif (scratch(istrt:istrt+5).eq.'common') then
              itype=-2
              istrt=istrt+6
            else
              itype=0
            endif
          endif
          if (itype.eq.-2) then
            do i=istrt,length
              if (scratch(i:i).eq.'/') scratch(i:i)=' '
            enddo
          endif
!
! Find number of words (other than reserved words above) and their
! starting and stopping indices; these are variable names. After
! that build global catalog. For each object discovered record its
! name, both the actual and small-case version for future match;
! length of the name, dimensions, if any (these are brackets
! immediately after the name), dimension length, and status. After
! that print output string. Note a line is printed out only if it
! is discovered that the object appear both in type and in common
! block declaration (hence varstat reaches 2), and it is the second
! appearance which triggers printing. Doing so immediately allows to
! preserve CPP-structure.
!
          nwords=0
          lswtch=.false.
          do i=istrt,length+1
            if (.not.lswtch .and. scratch(i:i).ge.'a'
     &                  .and. scratch(i:i).le.'z') then
              lswtch=.true.
              nwords=nwords+1
              is(nwords)=i
            elseif (lswtch .and. scratch(i:i).eq.' ') then
              lswtch=.false.
              ie(nwords)=i-1
            endif
          enddo
          is(nwords+1)=length+1
 
          do j=1,nwords
            io=0
            do i=is(j+1)-1,ie(j)+1,-1
              if (string(i:i).eq.'(') io=i
            enddo
            ic=0
            do i=ie(j)+1,is(j+1)-1,+1
              if (string(i:i).eq.')') ic=i
            enddo
            ivar=0
            do i=1,nvars
              if (ie(j)-is(j)+1 .eq. namelen(i)) then
                if (scratch(is(j):ie(j)).eq.
     &              svarname(i)(1:namelen(i))) ivar=i
              endif
            enddo
            if (ivar.eq.0) then
              ivar=nvars+1
              nvars=ivar
              varname(ivar) =string(is(j):ie(j))
              svarname(ivar)=scratch(is(j):ie(j))
              namelen(ivar)=ie(j)-is(j)+1
            endif
            if (io.gt.0 .and. ic.gt.0) then
              vardim(ivar)=string(io:ic)
              dimlen(ivar)=ic-io+1
            endif
            if (itype.gt.0 .or. itype.eq.-2) then
              varstat(ivar)=varstat(ivar)+1
            endif
 
 
            if (varstat(ivar).eq.2) then
              varstat(ivar)=0
c--#define ALLOCATE
#define ASSOCIATE
#if defined ALLOCATE
              if (dimlen(ivar).gt.0) then
                write(iout,'(9x,5A)') 'allocate(AMR_grids(level)%',
     &                            varname(ivar)(1:namelen(ivar)),
     &                            vardim(ivar)(1:dimlen(ivar)), ')'
              else
                write(iout,'(9x,5A)') 'allocate(AMR_grids(level)%',
     &                          varname(ivar)(1:namelen(ivar)), ')'
              endif
#elif defined ASSOCIATE
                write(iout,'(9x,5A)') varname(ivar)(1:namelen(ivar)),
     &       ' => AMR_grids(level)%', varname(ivar)(1:namelen(ivar))
#endif
            endif
          enddo
        elseif (string(istrt:istrt).eq.'#') then
          write(iout,'(A)') string(1:length)
        endif
      enddo
      stop
      end
 
