// Tests for tst.m. Run from the tests/ directory: magma -n test_all.m

load "../tst.m";

// --- TSTAssertEQ ---

procedure test_AssertEQ_integers()
    TSTAssertEQ(1, 1);
    TSTAssertEQ(0, 0);
    TSTAssertEQ(-5, -5);
end procedure;

procedure test_AssertEQ_strings()
    TSTAssertEQ("hello", "hello");
    TSTAssertEQ("", "");
end procedure;

procedure test_AssertEQ_sequences()
    TSTAssertEQ([1, 2, 3], [1, 2, 3]);
    TSTAssertEQ([], []);
end procedure;

procedure test_AssertEQ_fails_on_unequal()
    TSTAssertRaises(TSTAssertEQ, "", [1, 2]);
end procedure;

// --- TSTAssertNE ---

procedure test_AssertNE_integers()
    TSTAssertNE(1, 2);
    TSTAssertNE(0, -1);
end procedure;

procedure test_AssertNE_strings()
    TSTAssertNE("hello", "world");
end procedure;

procedure test_AssertNE_fails_on_equal()
    TSTAssertRaises(TSTAssertNE, "", [1, 1]);
end procedure;

// --- TSTAssertRaises ---

function raise_with_message(msg)
    error msg;
    return 0;
end function;

procedure test_AssertRaises_catches_expected_error()
    TSTAssertRaises(raise_with_message, "boom", ["boom"]);
end procedure;

procedure test_AssertRaises_matches_substring()
    TSTAssertRaises(raise_with_message, "sub", ["substring match"]);
end procedure;

function identity(x)
    return x;
end function;

procedure test_AssertRaises_fails_when_no_error_raised()
    TSTAssertRaises(TSTAssertRaises, "", [* identity, "", [42] *]);
end procedure;

procedure test_AssertRaises_fails_on_message_mismatch()
    TSTAssertRaises(TSTAssertRaises, "", [* raise_with_message, "expected", ["different"] *]);
end procedure;

procedure error_proc(x)
    error "proc boom";
end procedure;

procedure test_AssertRaises_with_procedure()
    TSTAssertRaises(error_proc, "proc boom", [42]);
end procedure;

// --- run all ---

test_AssertEQ_integers();
test_AssertEQ_strings();
test_AssertEQ_sequences();
test_AssertEQ_fails_on_unequal();
test_AssertNE_integers();
test_AssertNE_strings();
test_AssertNE_fails_on_equal();
test_AssertRaises_catches_expected_error();
test_AssertRaises_matches_substring();
test_AssertRaises_fails_when_no_error_raised();
test_AssertRaises_fails_on_message_mismatch();
test_AssertRaises_with_procedure();

print "All tests finished!";
exit;
