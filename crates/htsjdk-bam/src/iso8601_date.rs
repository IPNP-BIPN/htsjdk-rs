//! `Iso8601Date`, as `SAMTextHeaderCodec` uses it to rewrite a read group's `DT`.
//!
//! Ported from `htsjdk.samtools.util.DateParser.getCalendar` and `Iso8601Date`, htsjdk 5.0.0.
//! Parsing a `@RG` line replaces its `DT` with `decodeDate(dt).toString()`: the date is read as
//! ISO 8601 into a lenient UTC `GregorianCalendar`, truncated to the second, and formatted with
//! `yyyy-MM-dd'T'HH:mm:ssZ` in the JVM's time zone. The oracle runs in UTC, so
//! `2014-12-12T00:00:00-0500` becomes `2014-12-12T05:00:00+0000`.

/// `decodeDate(value).toString()` in a UTC JVM, or `None` where the value stays as written.
///
/// shortcut: `decodeDate`'s two fallbacks (`DateFormat.getDateTimeInstance().parse` and
/// `new Date(String)`) are not ported, so a non-ISO date stays as written; port them if a
/// header with such a date turns up. Years outside 1..=9999 also stay as written, since the
/// Julian calendar before 1582 and `yyyy`'s era handling are not reproduced either.
pub fn normalize(value: &str) -> Option<String> {
    let millis = parse(value)?;
    // `truncateMilliseconds` subtracts Java's `%`, which rounds a negative time toward zero.
    let seconds = (millis - millis % 1000) / 1000;
    let days = seconds.div_euclid(86_400);
    let second_of_day = seconds.rem_euclid(86_400);
    let (year, month, day) = civil_from_days(days);
    if !(1..=9999).contains(&year) {
        return None;
    }
    Some(format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}+0000",
        second_of_day / 3600,
        second_of_day / 60 % 60,
        second_of_day % 60
    ))
}

/// `DateParser.getCalendar(value).getTime()` in milliseconds, or `None` where it throws.
fn parse(value: &str) -> Option<i64> {
    let mut tokens = tokenize(value).into_iter().peekable();
    // `calendar.clear()` leaves every field at the epoch's.
    let (mut year, mut month, mut day) = (1970_i64, 0_i64, 1_i64);
    let (mut hour, mut minute, mut second, mut milli) = (0_i64, 0_i64, 0_i64, 0_i64);
    let mut offset_minutes = 0_i64;

    // `check`: false when nothing is left, an error when the next token is not the one asked.
    let check = |tokens: &mut std::iter::Peekable<std::vec::IntoIter<&str>>, want: &str| {
        match tokens.next() {
            None => Some(false),
            Some(t) if t == want => Some(true),
            Some(_) => None,
        }
    };

    'fields: {
        let Some(t) = tokens.next() else {
            break 'fields;
        };
        year = int(t)?;
        if !(check(&mut tokens, "-")? && tokens.peek().is_some()) {
            break 'fields;
        }
        month = int(tokens.next()?)? - 1;
        if !(check(&mut tokens, "-")? && tokens.peek().is_some()) {
            break 'fields;
        }
        day = int(tokens.next()?)?;
        if !(check(&mut tokens, "T")? && tokens.peek().is_some()) {
            break 'fields;
        }
        hour = int(tokens.next()?)?;
        if !(check(&mut tokens, ":")? && tokens.peek().is_some()) {
            break 'fields;
        }
        minute = int(tokens.next()?)?;
        let Some(mut tok) = tokens.next() else {
            break 'fields;
        };
        if tok == ":" {
            second = int(tokens.next()?)?;
            let Some(next) = tokens.next() else {
                break 'fields;
            };
            tok = next;
            if tok == "." {
                // The fraction is padded or cut to three digits, then read as milliseconds.
                let digits: String = format!("{:0<3}", tokens.next()?).chars().take(3).collect();
                milli = int(&digits)?;
                let Some(next) = tokens.next() else {
                    break 'fields;
                };
                tok = next;
            }
        }
        if tok != "Z" {
            let plus = match tok {
                "+" => true,
                "-" => false,
                _ => return None,
            };
            let mut tz_hour = int(tokens.next()?)?;
            let mut tz_minute = 0;
            if check(&mut tokens, ":")? && tokens.peek().is_some() {
                tz_minute = int(tokens.next()?)?;
            } else if tz_hour >= 100 {
                tz_minute = tz_hour % 100;
                tz_hour /= 100;
            }
            offset_minutes = if plus { 1 } else { -1 } * (tz_hour * 60 + tz_minute);
        }
    }

    // The calendar is lenient: an out-of-range field carries into the next one up.
    let months = year * 12 + month;
    let days = days_from_civil(months.div_euclid(12), months.rem_euclid(12) + 1) + day - 1;
    let seconds = days * 86_400 + hour * 3600 + (minute - offset_minutes) * 60 + second;
    Some(seconds * 1000 + milli)
}

/// `new StringTokenizer(value, "-T:.+Z", true)`: runs between delimiters, and each delimiter alone.
fn tokenize(value: &str) -> Vec<&str> {
    let mut tokens = Vec::new();
    let mut start = 0;
    for (i, c) in value.char_indices() {
        if "-T:.+Z".contains(c) {
            if start < i {
                tokens.push(&value[start..i]);
            }
            tokens.push(&value[i..i + 1]);
            start = i + 1;
        }
    }
    if start < value.len() {
        tokens.push(&value[start..]);
    }
    tokens
}

/// `Integer.parseInt` over ASCII digits, refused past `int`'s range.
fn int(token: &str) -> Option<i64> {
    if token.is_empty() || !token.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    token.parse::<i32>().ok().map(i64::from)
}

/// Days since 1970-01-01 of the first of `month` (1-based) in the proleptic Gregorian calendar.
fn days_from_civil(year: i64, month: i64) -> i64 {
    let y = if month <= 2 { year - 1 } else { year };
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let mp = (month + 9) % 12;
    let doy = (153 * mp + 2) / 5;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
}

/// The inverse of `days_from_civil`, with the day of the month.
fn civil_from_days(days: i64) -> (i64, i64, i64) {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    (year, month, day)
}

#[cfg(test)]
mod tests {
    use super::normalize;

    #[test]
    fn an_offset_date_is_rewritten_in_utc() {
        assert_eq!(
            normalize("2014-12-12T00:00:00-0500").as_deref(),
            Some("2014-12-12T05:00:00+0000")
        );
        assert_eq!(
            normalize("2014-12-12T00:00:00+05:30").as_deref(),
            Some("2014-12-11T18:30:00+0000")
        );
    }

    #[test]
    fn missing_fields_fall_back_to_the_epoch_and_fractions_are_dropped() {
        assert_eq!(
            normalize("2011-07-12").as_deref(),
            Some("2011-07-12T00:00:00+0000")
        );
        assert_eq!(
            normalize("2011").as_deref(),
            Some("2011-01-01T00:00:00+0000")
        );
        assert_eq!(
            normalize("2011-07-12T10:20:30.987Z").as_deref(),
            Some("2011-07-12T10:20:30+0000")
        );
        assert_eq!(
            normalize("2011-07-12T10:20Z").as_deref(),
            Some("2011-07-12T10:20:00+0000")
        );
    }

    #[test]
    fn the_calendar_is_lenient() {
        assert_eq!(
            normalize("2011-13-32").as_deref(),
            Some("2012-02-01T00:00:00+0000")
        );
    }

    #[test]
    fn an_unparseable_date_stays_as_written() {
        assert_eq!(normalize("Dec 12, 2014"), None);
        assert_eq!(normalize("2011/07/12"), None);
        assert_eq!(normalize("2011-07-12T10:20:30Q"), None);
    }
}
