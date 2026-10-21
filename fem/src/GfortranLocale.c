/*****************************************************************************
 *
 * TRAFOLO: random crashes in formatted output on MinGW.
 *
 * Without newlocale, libgfortran wraps every formatted transfer in
 *     old = setlocale(LC_NUMERIC, NULL); setlocale(LC_NUMERIC, "C");
 *     ... the transfer ...
 *     setlocale(LC_NUMERIC, old);
 * The UCRT rebuilds its locale data whenever setlocale sets a locale, also "C"
 * over "C", so `old` points to freed memory by the time it is restored. The
 * restore crashes once that memory has been released (random SIGSEGVs in long
 * runs, the more likely the more the heap is used), and when the freed bytes
 * start with a NUL it sets the user's locale instead of "C": C code between two
 * transfers (MATC, Lua) then parses and prints numbers with that locale's
 * decimal separator. elmer_gfortran_locale_fix points the setlocale import of
 * libgfortran at a wrapper that returns a copy of the queried name, valid until
 * the next query of that category. Only libgfortran's own calls go through it,
 * and the module holding the wrapper is pinned so that it stays loaded. It
 * needs the shared libgfortran (a static one does not link, see the dllimport
 * below). Elsewhere the function does nothing.
 *
 * Returns 1 when the workaround is in place, 0 where it is not needed (not
 * MinGW), 2 when no import of libgfortran resolves to setlocale and -1 when it
 * could not be installed; later calls return the first result. Call it first
 * from one thread; a concurrent first call patches the same slot to the same
 * wrapper.
 *
 *****************************************************************************/

#if defined(__MINGW32__)

#include <windows.h>
#include <locale.h>
#include <string.h>

typedef char *(__cdecl *SetlocaleFn)(int, const char *);

/* an entry point of libgfortran: its imported address lies in the libgfortran DLL */
__declspec(dllimport) void _gfortran_st_write(void *);

static SetlocaleFn real_setlocale = NULL;
static char locale_name[LC_MAX + 1][1024];

static char *__cdecl stable_setlocale(int category, const char *locale)
{
  char *name = real_setlocale(category, locale);

  if (locale != NULL || name == NULL || category < LC_MIN || category > LC_MAX)
    return name;
  strncpy(locale_name[category], name, sizeof(locale_name[category]) - 1);
  locale_name[category][sizeof(locale_name[category]) - 1] = '\0';
  return locale_name[category];
}

/* Point the import slot of module that resolves to target at stable_setlocale:
   1 when done (or done before), 0 when no slot resolves to target, -1 when the
   slot could not be written. Matching the resolved address, not the name,
   covers imports by name and by ordinal alike. */
static int patch_setlocale_import(HMODULE module, ULONG_PTR target)
{
  BYTE *base = (BYTE *)module;
  IMAGE_NT_HEADERS *nt = (IMAGE_NT_HEADERS *)(base + ((IMAGE_DOS_HEADER *)base)->e_lfanew);
  IMAGE_DATA_DIRECTORY dir = nt->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_IMPORT];
  IMAGE_IMPORT_DESCRIPTOR *imp;

  if (dir.VirtualAddress == 0) return 0;
  for (imp = (IMAGE_IMPORT_DESCRIPTOR *)(base + dir.VirtualAddress); imp->Name != 0; imp++) {
    IMAGE_THUNK_DATA *slot;

    for (slot = (IMAGE_THUNK_DATA *)(base + imp->FirstThunk); slot->u1.Function != 0; slot++) {
      DWORD protect;

      if (slot->u1.Function == (ULONG_PTR)stable_setlocale) return 1;
      if (slot->u1.Function != target) continue;
      /* the import table is read-only */
      if (!VirtualProtect(&slot->u1.Function, sizeof(slot->u1.Function), PAGE_READWRITE, &protect))
        return -1;
      real_setlocale = (SetlocaleFn)target;
      slot->u1.Function = (ULONG_PTR)stable_setlocale;
      VirtualProtect(&slot->u1.Function, sizeof(slot->u1.Function), protect, &protect);
      return 1;
    }
  }
  return 0;
}

int elmer_gfortran_locale_fix(void)
{
  static int tried = 0, status = -1;
  HMODULE gfortran, self;

  if (tried) return status;
  tried = 1;
  /* libgfortran itself, not a module that merely calls it, and this module
     pinned: libgfortran's import will point into it */
  if (GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                         GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                         (LPCSTR)(ULONG_PTR)&_gfortran_st_write, &gfortran) &&
      GetProcAddress(gfortran, "_gfortran_st_write") != NULL &&
      GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_PIN,
                         (LPCSTR)(ULONG_PTR)&stable_setlocale, &self)) {
    /* the C runtime's setlocale, as this module resolves it */
    status = patch_setlocale_import(gfortran, (ULONG_PTR)&setlocale);
    if (status == 0) status = 2;
  }
  return status;
}

#else

int elmer_gfortran_locale_fix(void)
{
  return 0;
}

#endif
