# magma-unittest

A minimal unit testing framework for [Magma](http://magma.maths.usyd.edu.au/magma/).

## Usage

Add this repository as a git submodule, then load `tst.m` from your test entry point:

    load "path/to/magma-unittest/tst.m";

See `tests/test_all.m` for working examples of every function in the API.

## API

- `TSTAssertEQ(computed, expected)` -- asserts `computed eq expected`, aborts on failure
- `TSTAssertNE(computed, expected)` -- asserts `computed ne expected`, aborts on failure
- `TSTAssertRaises(f, message, args)` -- asserts that calling `f` with argument list `args`
  raises an error whose message contains the substring `message`

## Adding as a submodule

    git submodule add git@github.com:nt-lib/magma-unittest.git tests/magma-unittest

Then in your test entry point:

    load "magma-unittest/tst.m";

## Testing

To run the tests for this repository, you need a working Magma installation.

If Magma is available locally:

    make test

If Magma is only available on a remote host over SSH:

    make test_remote ssh=<hostname>

This rsyncs the repository to `/tmp/magma-unittest` on the remote host and runs
the tests there.

## Licence

Copyright (C) 2025 Maarten Derickx.
Licensed under the GNU General Public Licence v3. See LICENSE.
