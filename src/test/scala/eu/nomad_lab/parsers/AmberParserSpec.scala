/*
   Copyright 2016-2017 The NOMAD Developers Group

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
 */
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
