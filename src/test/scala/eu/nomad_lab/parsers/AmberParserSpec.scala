package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object AmberParserSpec extends Specification {
  "AmberParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(AmberParser, "parsers/amber/test/examples/03_Prod.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(AmberParser, "parsers/amber/test/examples/03_Prod.out", "json") must_== ParseResult.ParseSuccess
    }
  }
}
